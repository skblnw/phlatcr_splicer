#!/usr/bin/env python3
"""
Unified TCR-pMHC Structure Analyzer

A unified analyzer that handles both MHC-I and MHC-II complexes in a single system.
This analyzer automatically identifies and classifies all chain types, performs
spatial clustering to group complexes, and pairs chains based on proximity.

Key Features:
- Handles both MHC-I (with B2M) and MHC-II (α/β heterodimer) complexes
- Uses sequence alignment for high-confidence TCR α/β discrimination
- Performs DBSCAN clustering for multi-complex detection
- Supports both PDB and CIF file formats
- Proximity-based reclassification of unknown chains

Author: skblnw
"""

import os
import logging
import warnings
from typing import Dict, List, Tuple, Optional, Set
from collections import defaultdict
import numpy as np
import re

try:
    from Bio.PDB import PDBParser, MMCIFParser, Structure, Model, Chain
    from Bio.PDB.MMCIF2Dict import MMCIF2Dict
    from Bio.SeqUtils import seq1
    from Bio.Seq import Seq
    from Bio import Align
except ImportError:
    raise ImportError("BioPython is required. Install with: pip install biopython")

try:
    from sklearn.cluster import DBSCAN
except ImportError:
    raise ImportError("scikit-learn is required. Install with: pip install scikit-learn")

# Configure logging
logger = logging.getLogger(__name__)

# Reference sequences for TCR constant regions (conserved domains)
TCR_ALPHA_CONST_HUMAN = "IQNPDPAVYQLRDSKSSDKSVCLFTDFDSQTNVSQSKDSDVYITDKTVLDMRSMDFKSNSAVAWSNKSDFACANAFNNSIIPEDTFFPSPESSCDVKLVEKSFETDTNLNFQNLSVIGFRILLLKVAGFNLLMTLRLWSS"
TCR_BETA_CONST_HUMAN = "CASSLGGNQPQHFGDGTRLSILEDLNKVFPPEVAVFEPSEAEISHTQKATLVCLATGFFPDHVELSWWVNGKEVHSGVSTDPQPLKEQPALNDSRYCLSSRLRVSATFWQNPRNHFRCQVQFYGLSENDEWTQDRAKPVTQIVSAEAWGRADCGFTSVSYQQGVLSATILYEILLGKATLYAVLVSALVLMAMVKRKDF"
TCR_ALPHA_CONST_MOUSE = "MAMLELLTLAFLGIWAFNQRADGQIYNQEPATNENISEATYNASLCSLTESGKSYFFWYKQEPGAGLQLLTYIFSNMDMKQDQRLTVLLLNKKDKHLSLRIADTQTGDSAIYFCAVSGGYQKVTFGIGTKLQVIP"
TCR_BETA_CONST_MOUSE = "MGTQILCLMLIFLGSGNTGAGVTQTPKFQVLKTGQSMTLQCAQDMNHEYMSWYRQDPGMGLRLIHYSVGAGITDQGEVPNGYNVSRSTTEDFPLRLLSAAPSQTSVYFCASSYVGNTGELFFGEGTRLTVV"

# Unified chain type configuration supporting both MHC-I and MHC-II
CHAIN_TYPE_CONFIG = {
    'TCR_ALPHA': {
        'length_range': (90, 320),
        'patterns': [
            r'CA[A-Z]{1,25}FG[A-Z]{2}GT',  # General CDR3+J for alpha
            r'DSAGY',                       # Conserved motif
            r'L[IV]LLKVAGFNLLMTLRLWSS',     # Transmembrane/cytoplasmic
            r'ANAFNNSIIPEDTFFPSPESS',       # Constant region motif
            r'Q[A-Z]G[A-Z]{2}AT',           # Pattern for broader coverage
            r'DLYSGLPS',                    # Additional TCR-α motif
            r'FGNGTQL',                     # J-region pattern
        ],
        'weight': 2.5,
        'ref_seqs': [TCR_ALPHA_CONST_HUMAN, TCR_ALPHA_CONST_MOUSE]
    },
    'TCR_BETA': {
        'length_range': (90, 320),
        'patterns': [
            r'CASS[A-Z]{1,25}FG[A-Z]{2}GT', # General CDR3+J for beta
            r'DLNKVFPPEVAVFE',              # Constant region
            r'CLATGFFPDHVELSWWVNGKEVHS',    # Constant region
            r'ATFWQNPRNHFRCQVQFYGL',        # Constant region
            r'EDLKNVFPPEVAVFE',              # Alternative constant region
            r'FGPGTRL',                      # J-region pattern
        ],
        'weight': 2.5,
        'ref_seqs': [TCR_BETA_CONST_HUMAN, TCR_BETA_CONST_MOUSE]
    },
    'MHC_I_ALPHA': {
        'length_range': (240, 300),
        'patterns': [
            r'GSHSMR',      # MHC-I specific N-terminal
            r'YFYT',        # Alpha1 domain signature
            r'HLA-[ABC]',   # HLA class I
            r'WMEYE',       # Common MHC-I motif
            r'QEYR',        # MHC-I pattern
            r'EPRAP',       # N-terminal region
            r'MRYFYT',      # Extended pattern
        ],
        'weight': 1.0
    },
    'B2M': {
        'length_range': (80, 120),
        'patterns': [
            r'IQRTPKIQ',    # B2M specific
            r'B2M',         # Direct mention
            r'VNHVTLS',     # B2M conserved
            r'KSNFLNC',     # B2M pattern
        ],
        'weight': 1.0
    },
    'MHC_II_ALPHA': {
        'length_range': (150, 220),
        'patterns': [
            r'WRLEEF',      # MHC-II α specific
            r'LREPNV',      # MHC-II α pattern
            r'WLRNGK',      # MHC-II α motif
            r'FLPRED',      # MHC-II α signature
            r'[RK]G[RK]G',  # MHC-specific motif
            r'GADGV',       # MHC-II alpha pattern
            r'[VL]P[RK]',   # Additional pattern
            r'MGSGWV',      # N-terminal signature
        ],
        'weight': 1.0
    },
    'MHC_II_BETA': {
        'length_range': (150, 230),
        'patterns': [
            r'NGTER',       # MHC-II β specific
            r'FDSDVG',      # MHC-II β pattern
            r'YWNSQK',      # MHC-II β motif
            r'VDTYCR',      # MHC-II β signature
            r'GFYPGS',      # MHC-II β pattern
            r'[VL]G[DE]',   # Additional pattern
            r'Y[FY]',       # Short pattern
            r'RATPE',       # β1 domain pattern
        ],
        'weight': 1.0
    },
    'PEPTIDE': {
        'length_range': (5, 40),
        'patterns': [
            r'^[ACDEFGHIKLMNPQRSTVWY]{5,40}$',  # Full peptide sequence
        ],
        'weight': 2.0
    },
}


class TCRpMHCAnalyzer:
    """
    Unified analyzer for TCR-pMHC complex structures.
    
    Handles both MHC-I and MHC-II complexes, automatically identifying chain types,
    clustering complexes, and pairing chains based on spatial proximity.
    """
    
    def __init__(self, verbose: bool = False, eps: float = 60.0,
                 align_score: int = 20, align_ratio: float = 1.5,
                 tcr_pair_dist: float = 50.0, mhc1_pair_dist: float = 50.0, 
                 mhc2_pair_dist: float = 50.0, pep_mhc1_dist: float = 40.0, 
                 pep_mhc2_dist: float = 60.0, tcr_pmhc_dist: float = 150.0,
                 tcr_pmhc_dist_unpaired: float = 120.0, reclassify_dist: float = 50.0):
        """
        Initialize the unified analyzer with configurable parameters.
        
        Args:
            verbose: Enable verbose logging for debugging
            eps: DBSCAN clustering distance in Angstroms
            align_score: Minimum alignment score for confident TCR classification
            align_ratio: Score ratio to distinguish TCR-alpha vs TCR-beta
            tcr_pair_dist: Max distance for TCR α/β pairing
            mhc1_pair_dist: Max distance for MHC-I/B2M pairing
            mhc2_pair_dist: Max distance for MHC-II α/β pairing
            pep_mhc1_dist: Max distance for Peptide/MHC-I pairing
            pep_mhc2_dist: Max distance for Peptide/MHC-II pairing
            tcr_pmhc_dist: Max distance for TCR/pMHC pairing
            tcr_pmhc_dist_unpaired: Max distance for single-chain TCR/pMHC pairing
            reclassify_dist: Max distance to reclassify UNKNOWN chains by proximity
        """
        self.verbose = verbose
        self.eps = eps
        self.align_score = align_score
        self.align_ratio = align_ratio
        self.tcr_pair_dist = tcr_pair_dist
        self.mhc1_pair_dist = mhc1_pair_dist
        self.mhc2_pair_dist = mhc2_pair_dist
        self.pep_mhc1_dist = pep_mhc1_dist
        self.pep_mhc2_dist = pep_mhc2_dist
        self.tcr_pmhc_dist = tcr_pmhc_dist
        self.tcr_pmhc_dist_unpaired = tcr_pmhc_dist_unpaired
        self.reclassify_dist = reclassify_dist
        
        # Initialize parsers and aligner
        self.pdb_parser = PDBParser(QUIET=True)
        self.cif_parser = MMCIFParser(QUIET=True)
        self.aligner = Align.PairwiseAligner()
        self.aligner.mode = 'local'
        self.aligner.match_score = 2
        self.aligner.mismatch_score = -1
        self.aligner.open_gap_score = -5
        self.aligner.extend_gap_score = -0.5
        
        if verbose:
            logger.setLevel(logging.DEBUG)
        else:
            logger.setLevel(logging.INFO)
    
    def _get_chain_sequence(self, chain) -> str:
        """Extract the amino acid sequence from a PDB chain object."""
        seq = []
        for res in chain:
            # Skip non-standard residues (water, ligands)
            if res.id[0] != ' ':
                continue
            try:
                # Convert 3-letter to 1-letter amino acid code
                seq.append(seq1(res.resname))
            except KeyError:
                # Use 'X' for unknown residues
                seq.append('X')
        return ''.join(seq)
    
    def _identify_chain_type(self, seq: str) -> str:
        """
        Identify chain type using a three-tier hierarchical approach:
        1. Check for peptides by length
        2. Attempt high-confidence TCR classification by sequence alignment
        3. Fall back to comprehensive pattern matching for all types
        """
        length = len(seq)
        logger.debug(f"Identifying chain (len={length}): {seq[:50]}...")
        
        # Priority 1: Classify peptides based on length
        if CHAIN_TYPE_CONFIG['PEPTIDE']['length_range'][0] <= length <= CHAIN_TYPE_CONFIG['PEPTIDE']['length_range'][1]:
            logger.debug("Classified as PEPTIDE based on length")
            return 'PEPTIDE'
        
        # Priority 2: Attempt high-confidence TCR classification via alignment
        is_potential_tcr = (CHAIN_TYPE_CONFIG['TCR_ALPHA']['length_range'][0] <= length <= 
                           CHAIN_TYPE_CONFIG['TCR_ALPHA']['length_range'][1])
        
        if is_potential_tcr:
            logger.debug("Potential TCR chain, attempting alignment")
            
            # Align against reference TCR constant region sequences
            alpha_scores = []
            beta_scores = []
            
            for ref in CHAIN_TYPE_CONFIG['TCR_ALPHA']['ref_seqs']:
                alignments = self.aligner.align(seq, ref)
                if alignments:
                    alpha_scores.append(alignments[0].score)
            
            for ref in CHAIN_TYPE_CONFIG['TCR_BETA']['ref_seqs']:
                alignments = self.aligner.align(seq, ref)
                if alignments:
                    beta_scores.append(alignments[0].score)
            
            alpha_score = max(alpha_scores) if alpha_scores else 0
            beta_score = max(beta_scores) if beta_scores else 0
            logger.debug(f"Alignment scores -> Alpha: {alpha_score}, Beta: {beta_score}")
            
            # Classify if one score is significantly higher
            if alpha_score > self.align_score and alpha_score > beta_score * self.align_ratio:
                return 'TCR_ALPHA'
            if beta_score > self.align_score and beta_score > alpha_score * self.align_ratio:
                return 'TCR_BETA'
            
            # Handle ambiguous cases where scores are close
            if (alpha_score > self.align_score and beta_score > self.align_score and 
                abs(alpha_score - beta_score) / max(alpha_score, beta_score) < 0.1):
                return 'TCR_ALPHA'  # Default to Alpha as tie-breaker
            
            logger.debug("Alignment inconclusive, falling back to pattern matching")
        
        # Priority 3: Use pattern matching for all chain types
        scores = defaultdict(float)
        for chain_type, config in CHAIN_TYPE_CONFIG.items():
            if chain_type == 'PEPTIDE':
                continue
            
            lmin, lmax = config['length_range']
            if not (lmin <= length <= lmax):
                continue
            
            pattern_score = sum(1 for pattern in config['patterns'] if re.search(pattern, seq))
            scores[chain_type] = pattern_score * config.get('weight', 1.0)
        
        logger.debug(f"Pattern scores: {dict(scores)}")
        
        if not scores or max(scores.values()) == 0:
            logger.debug("No candidate found, classified as UNKNOWN")
            return 'UNKNOWN'
        
        # Return the best candidate from pattern matching
        best_candidate = max(scores, key=scores.get)
        logger.debug(f"Classified as {best_candidate} with score {scores[best_candidate]}")
        return best_candidate
    
    def _compute_chain_center(self, chain) -> np.ndarray:
        """Compute the geometric center (centroid) of a chain using C-alpha atoms."""
        coords = [atom.coord for res in chain for atom in res if atom.name == 'CA']
        if not coords:
            return np.array([np.nan, np.nan, np.nan])
        return np.mean(coords, axis=0)
    
    def analyze_pdb(self, structure_file: str) -> Tuple[Dict[str, Dict], Dict[str, str]]:
        """
        Analyze a PDB or CIF file to identify chain types and complexes.
        
        Args:
            structure_file: Path to PDB or CIF file
            
        Returns:
            Tuple of (results_dict, chain_id_map)
            - results_dict: Dictionary of complexes with chain assignments and pairings
            - chain_id_map: Mapping from internal to author chain IDs (for CIF files)
        """
        chain_id_map = {}
        
        # Handle CIF files with chain ID mapping
        if structure_file.lower().endswith('.cif'):
            try:
                mmcif_dict = MMCIF2Dict(structure_file)
                label_ids = mmcif_dict.get('_struct_asym.label_asym_id', [])
                auth_ids = mmcif_dict.get('_struct_asym.auth_asym_id', [])
                
                if isinstance(label_ids, str):
                    label_ids = [label_ids]
                if isinstance(auth_ids, str):
                    auth_ids = [auth_ids]
                
                if label_ids and auth_ids and len(label_ids) == len(auth_ids):
                    chain_id_map = dict(zip(label_ids, auth_ids))
            except (KeyError, ValueError):
                # If CIF parsing fails, continue without chain mapping
                pass
        
        # Parse structure
        parser = self.cif_parser if structure_file.lower().endswith('.cif') else self.pdb_parser
        structure = parser.get_structure('struct', structure_file)
        
        # Check if structure has models
        if len(structure) == 0:
            logger.warning(f"No models found in {structure_file}")
            return {}, chain_id_map
        
        chain_centers = {}
        chain_types = {}
        chain_ids = []
        
        # Process all chains
        model = structure[0]
        for chain in model:
            cid = chain.id
            seq = self._get_chain_sequence(chain)
            
            if len(seq) < 5:
                continue
            
            chain_types[cid] = self._identify_chain_type(seq)
            
            if chain_types[cid] == 'UNKNOWN' and self.verbose:
                logger.debug(f"Unknown chain {cid} sequence: {seq}")
            
            chain_centers[cid] = self._compute_chain_center(chain)
            chain_ids.append(cid)
        
        if not chain_centers:
            logger.warning(f"No valid chains found in {structure_file}")
            return {}, chain_id_map
        
        # Spatial clustering using DBSCAN
        centers = np.array([chain_centers[cid] for cid in chain_ids if not np.isnan(chain_centers[cid]).any()])
        valid_chain_ids = [cid for cid in chain_ids if not np.isnan(chain_centers[cid]).any()]
        
        if len(centers) < 2:
            clusters = np.zeros(len(centers), dtype=int) if len(centers) > 0 else []
        else:
            db = DBSCAN(eps=self.eps, min_samples=1, metric='euclidean')
            clusters = db.fit_predict(centers)
        
        # Group chains by cluster and ensure MHC-II chains are paired
        complex_groups = defaultdict(list)
        for idx, label in enumerate(clusters):
            chain_id = valid_chain_ids[idx]
            complex_groups[label].append(chain_id)
        
        # Ensure MHC-II alpha and beta chains are in the same complex
        for label, chains in list(complex_groups.items()):
            mhc_alpha_in_group = any(chain_types[c] == 'MHC_II_ALPHA' for c in chains)
            mhc_beta_in_group = any(chain_types[c] == 'MHC_II_BETA' for c in chains)
            
            # If only one type of MHC-II chain is present, find its pair
            if mhc_alpha_in_group ^ mhc_beta_in_group:  # XOR: only one is present
                for other_label, other_chains in complex_groups.items():
                    if label == other_label:
                        continue
                    
                    if mhc_alpha_in_group and any(chain_types[c] == 'MHC_II_BETA' for c in other_chains):
                        # Move MHC-II beta chain to this group
                        beta_chains = [c for c in other_chains if chain_types[c] == 'MHC_II_BETA']
                        complex_groups[label].extend(beta_chains)
                        complex_groups[other_label] = [c for c in other_chains if chain_types[c] != 'MHC_II_BETA']
                        if not complex_groups[other_label]:
                            del complex_groups[other_label]
                        break
                    
                    elif mhc_beta_in_group and any(chain_types[c] == 'MHC_II_ALPHA' for c in other_chains):
                        # Move MHC-II alpha chain to this group
                        alpha_chains = [c for c in other_chains if chain_types[c] == 'MHC_II_ALPHA']
                        complex_groups[label].extend(alpha_chains)
                        complex_groups[other_label] = [c for c in other_chains if chain_types[c] != 'MHC_II_ALPHA']
                        if not complex_groups[other_label]:
                            del complex_groups[other_label]
                        break
        
        # Reclassify UNKNOWN chains based on proximity
        self._reclassify_unknown_chains(chain_types, chain_centers, structure)
        
        # Pair chains within each group
        results = {}
        for group_id, chains in complex_groups.items():
            group_types = {cid: chain_types[cid] for cid in chains}
            pairs = self._pair_chains_in_group(structure, chains, group_types, chain_centers)
            results[f'complex_{group_id + 1}'] = {
                'chains': group_types,
                'pairs': pairs
            }
        
        return results, chain_id_map
    
    def _reclassify_unknown_chains(self, chain_types: Dict[str, str], 
                                   chain_centers: Dict[str, np.ndarray],
                                   structure: Structure):
        """Reclassify UNKNOWN chains based on proximity to known chains."""
        
        # Get chains of each type
        b2m_chains = [k for k, v in chain_types.items() if v == 'B2M']
        mhc_beta_chains = [k for k, v in chain_types.items() if v == 'MHC_II_BETA']
        mhc_alpha_chains = [k for k, v in chain_types.items() if v == 'MHC_II_ALPHA']
        mhc1_alpha_chains = [k for k, v in chain_types.items() if v == 'MHC_I_ALPHA']
        
        for cid in list(chain_types.keys()):
            if chain_types[cid] == 'UNKNOWN':
                # Rule 1: UNKNOWN near MHC-I alpha is likely B2M
                for mhc1a_cid in mhc1_alpha_chains:
                    if np.linalg.norm(chain_centers[cid] - chain_centers[mhc1a_cid]) < self.reclassify_dist:
                        chain_types[cid] = 'B2M'
                        break
                
                # Rule 2: UNKNOWN near B2M is likely MHC-I alpha
                if chain_types[cid] == 'UNKNOWN':
                    for b2m_cid in b2m_chains:
                        if np.linalg.norm(chain_centers[cid] - chain_centers[b2m_cid]) < self.reclassify_dist:
                            chain_types[cid] = 'MHC_I_ALPHA'
                            break
                
                # Rule 3: UNKNOWN near MHC-II beta is likely MHC-II alpha
                if chain_types[cid] == 'UNKNOWN':
                    for mhc2b_cid in mhc_beta_chains:
                        if np.linalg.norm(chain_centers[cid] - chain_centers[mhc2b_cid]) < self.reclassify_dist:
                            chain_types[cid] = 'MHC_II_ALPHA'
                            break
                
                # Rule 4: UNKNOWN near MHC-II alpha is likely MHC-II beta
                if chain_types[cid] == 'UNKNOWN':
                    for mhc2a_cid in mhc_alpha_chains:
                        if np.linalg.norm(chain_centers[cid] - chain_centers[mhc2a_cid]) < self.reclassify_dist:
                            chain_types[cid] = 'MHC_II_BETA'
                            break
                
                # Rule 5: Fallback for UNKNOWN chains in MHC length range
                if chain_types[cid] == 'UNKNOWN':
                    length = len(self._get_chain_sequence(structure[0][cid]))
                    if 150 <= length <= 230:
                        chain_types[cid] = 'MHC_II_ALPHA'
    
    def _pair_chains_in_group(self, structure: Structure, chains: List[str], 
                              types: Dict[str, str], centers: Dict[str, np.ndarray]) -> List[Dict]:
        """
        Group chains into TCR-pMHC complexes based on proximity.
        
        Returns a list of paired complexes with their component chains.
        """
        # Extract chains of each type
        tcr_alpha = [c for c in chains if types[c] == 'TCR_ALPHA']
        tcr_beta = [c for c in chains if types[c] == 'TCR_BETA']
        mhc1_alpha = [c for c in chains if types[c] == 'MHC_I_ALPHA']
        b2m = [c for c in chains if types[c] == 'B2M']
        mhc2_alpha = [c for c in chains if types[c] == 'MHC_II_ALPHA']
        mhc2_beta = [c for c in chains if types[c] == 'MHC_II_BETA']
        peptides = [c for c in chains if types[c] == 'PEPTIDE']
        
        # 1. Form TCR heterodimers (αβ)
        tcr_dimers = []
        used_alpha, used_beta = set(), set()
        
        for alpha in tcr_alpha:
            if tcr_beta:
                distances = [(b, np.linalg.norm(centers[alpha] - centers[b])) 
                            for b in tcr_beta if b not in used_beta]
                if distances:
                    closest_beta, min_dist = min(distances, key=lambda x: x[1])
                    if min_dist < self.tcr_pair_dist:
                        tcr_dimers.append({
                            'alpha': alpha, 
                            'beta': closest_beta,
                            'center': np.mean([centers[alpha], centers[closest_beta]], axis=0)
                        })
                        used_alpha.add(alpha)
                        used_beta.add(closest_beta)
        
        # Add unpaired TCR chains
        for alpha in tcr_alpha:
            if alpha not in used_alpha:
                tcr_dimers.append({'alpha': alpha, 'beta': None, 'center': centers[alpha]})
        for beta in tcr_beta:
            if beta not in used_beta:
                tcr_dimers.append({'alpha': None, 'beta': beta, 'center': centers[beta]})
        
        # 2. Form pMHC complexes
        pmhc_complexes = []
        used_peptides = set()
        
        # MHC-I complexes (MHC-I-alpha + B2M + optional Peptide)
        used_b2m = set()
        for m1a in mhc1_alpha:
            if b2m:
                distances = [(b, np.linalg.norm(centers[m1a] - centers[b]))
                            for b in b2m if b not in used_b2m]
                if distances:
                    closest_b2m, min_dist_b2m = min(distances, key=lambda x: x[1])
                    if min_dist_b2m < self.mhc1_pair_dist:
                        mhc_center = np.mean([centers[m1a], centers[closest_b2m]], axis=0)
                        
                        # Find closest peptide (optional)
                        peptide_chain = None
                        complex_center = mhc_center
                        
                        if peptides:
                            pep_distances = [(p, np.linalg.norm(mhc_center - centers[p]))
                                           for p in peptides if p not in used_peptides]
                            if pep_distances:
                                closest_pep, min_dist_pep = min(pep_distances, key=lambda x: x[1])
                                if min_dist_pep < self.pep_mhc1_dist:
                                    peptide_chain = closest_pep
                                    complex_center = np.mean([mhc_center, centers[peptide_chain]], axis=0)
                                    used_peptides.add(peptide_chain)
                        
                        pmhc_complexes.append({
                            'type': 'MHC-I',
                            'mhc_alpha': m1a,
                            'b2m': closest_b2m,
                            'peptide': peptide_chain,
                            'center': complex_center
                        })
                        used_b2m.add(closest_b2m)
        
        # MHC-II complexes (MHC-II-alpha + MHC-II-beta + optional Peptide)
        used_mhc2_beta = set()
        for m2a in mhc2_alpha:
            if mhc2_beta:
                distances = [(b, np.linalg.norm(centers[m2a] - centers[b]))
                            for b in mhc2_beta if b not in used_mhc2_beta]
                if distances:
                    closest_m2b, min_dist_m2b = min(distances, key=lambda x: x[1])
                    if min_dist_m2b < self.mhc2_pair_dist:
                        mhc_center = np.mean([centers[m2a], centers[closest_m2b]], axis=0)
                        
                        # Find closest peptide (optional)
                        peptide_chain = None
                        complex_center = mhc_center
                        
                        if peptides:
                            pep_distances = [(p, np.linalg.norm(mhc_center - centers[p]))
                                           for p in peptides if p not in used_peptides]
                            if pep_distances:
                                closest_pep, min_dist_pep = min(pep_distances, key=lambda x: x[1])
                                if min_dist_pep < self.pep_mhc2_dist:
                                    peptide_chain = closest_pep
                                    complex_center = np.mean([mhc_center, centers[peptide_chain]], axis=0)
                                    used_peptides.add(peptide_chain)
                        
                        pmhc_complexes.append({
                            'type': 'MHC-II',
                            'mhc_alpha': m2a,
                            'mhc_beta': closest_m2b,
                            'peptide': peptide_chain,
                            'center': complex_center
                        })
                        used_mhc2_beta.add(closest_m2b)
        
        # 3. Pair TCR dimers with pMHC complexes
        final_complexes = []
        used_pmhc_indices, used_tcr_indices = set(), set()
        
        # Prioritize complete TCR dimers (αβ)
        for i, tcr in enumerate(tcr_dimers):
            if tcr.get('alpha') and tcr.get('beta'):
                if pmhc_complexes:
                    distances = [(pmhc, np.linalg.norm(tcr['center'] - pmhc['center']), j)
                                for j, pmhc in enumerate(pmhc_complexes) 
                                if j not in used_pmhc_indices]
                    if distances:
                        closest_pmhc, min_dist, closest_idx = min(distances, key=lambda x: x[1])
                        if min_dist < self.tcr_pmhc_dist:
                            final_complexes.append({
                                'tcr_alpha': tcr.get('alpha'),
                                'tcr_beta': tcr.get('beta'),
                                **closest_pmhc
                            })
                            used_pmhc_indices.add(closest_idx)
                            used_tcr_indices.add(i)
        
        # Pair single TCR chains with remaining pMHCs
        for i, tcr in enumerate(tcr_dimers):
            if i in used_tcr_indices:
                continue
            if tcr.get('alpha') or tcr.get('beta'):
                if pmhc_complexes:
                    distances = [(pmhc, np.linalg.norm(tcr['center'] - pmhc['center']), j)
                                for j, pmhc in enumerate(pmhc_complexes)
                                if j not in used_pmhc_indices]
                    if distances:
                        closest_pmhc, min_dist, closest_idx = min(distances, key=lambda x: x[1])
                        if min_dist < self.tcr_pmhc_dist_unpaired:
                            final_complexes.append({
                                'tcr_alpha': tcr.get('alpha'),
                                'tcr_beta': tcr.get('beta'),
                                **closest_pmhc
                            })
                            used_pmhc_indices.add(closest_idx)
                            used_tcr_indices.add(i)
        
        # Add orphan pMHCs
        for i, pmhc in enumerate(pmhc_complexes):
            if i not in used_pmhc_indices:
                final_complexes.append({
                    'tcr_alpha': None,
                    'tcr_beta': None,
                    **pmhc
                })
        
        # Add orphan TCRs
        for i, tcr in enumerate(tcr_dimers):
            if i not in used_tcr_indices:
                final_complexes.append({
                    'tcr_alpha': tcr.get('alpha'),
                    'tcr_beta': tcr.get('beta'),
                    'type': 'UNPAIRED_TCR'
                })
        
        return final_complexes