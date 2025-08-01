"""
MHC-II Complex Analyzer

A specialized tool for identifying and classifying chains in MHC-II-TCR complexes.
Unlike MHC-I complexes, MHC-II complexes contain:
- MHC-II α chain (~180-200 residues)  
- MHC-II β chain (~190-210 residues)
- Longer peptide (12-25 residues)
- TCR α chain (180-230 residues)
- TCR β chain (230-290 residues)
- No β2-microglobulin

Author: AI Assistant
"""

from typing import Dict, List, Tuple, Optional, Set
from collections import defaultdict
import numpy as np

try:
    from Bio.PDB import PDBParser, Structure, Structure, Model, Chain, Residue
    from Bio.PDB.PDBIO import PDBIO
    from Bio.SeqUtils import ProtParam
    from Bio.Seq import Seq
except ImportError:
    print("Error: BioPython is required. Install with: pip install biopython")
    raise

import os
import tempfile


class pMHCIITCRAnalyzer:
    """
    Analyzer for MHC-II-TCR complex structures.
    
    Identifies and classifies the following chain types:
    - mhc_ii_alpha: MHC-II α chain
    - mhc_ii_beta: MHC-II β chain  
    - peptide: Presented peptide antigen (longer than MHC-I)
    - tcr_alpha: T-cell receptor α chain
    - tcr_beta: T-cell receptor β chain
    """
    
    def __init__(self, verbose: bool = False):
        """
        Initialize the MHC-II analyzer.
        
        Args:
            verbose: If True, print detailed analysis information
        """
        self.verbose = verbose
        self.parser = PDBParser(QUIET=True)
        
        # Load patterns for each chain type
        self.mhc_ii_patterns = self._load_mhc_ii_patterns()
        self.tcr_patterns = self._load_tcr_patterns()
        self.peptide_patterns = self._load_peptide_patterns()
        
    def _load_mhc_ii_patterns(self) -> Dict:
        """Load MHC-II alpha and beta chain patterns and characteristics."""
        return {
            'alpha': {
                'length_range': (170, 210),  # MHC-II α chains
                'optimal_length': (180, 200),
                'highly_specific_motifs': [
                    'MGSGWV',    # Common MHC-II α N-terminal
                    'MILNKA',    # MHC-II α leader sequence variant
                    'RFDSDK',    # α1 domain pattern
                    'EPRAPW',    # α1 domain signature
                    'QFRVD',     # α1-α2 linker region
                ],
                'conserved_motifs': [
                    'HIVQR',     # α1 domain
                    'WVLVL',     # α2 domain
                    'HSLTE',     # α2 domain
                    'VTQQP',     # α2 domain
                    'CRHNY',     # α3 domain
                ],
                'n_terminal_signatures': [
                    'MGSGWV',    # Very specific α chain start
                    'MILNKA',    # Alternative start
                    'GSGWV',     # Variant without Met
                ],
                'c_terminal_patterns': [
                    'CRHNY',     # α3 domain
                    'HVTLS',     # C-terminal region
                    'GACIF',     # Transmembrane region
                ],
                'domain_boundaries': {
                    'alpha1': (1, 85),
                    'alpha2': (86, 175),
                    'alpha3': (176, 200),
                },
            },
            'beta': {
                'length_range': (180, 220),  # MHC-II β chains  
                'optimal_length': (190, 210),
                'highly_specific_motifs': [
                    'MSWKKA',    # Common MHC-II β N-terminal
                    'MIHSPS',    # MHC-II β leader sequence variant
                    'PVTVY',     # β1 domain signature
                    'NYGVV',     # β1 domain pattern
                    'WNGLR',     # β1-β2 linker
                ],
                'conserved_motifs': [
                    'FLNGQ',     # β1 domain
                    'YVGC',      # β1 domain
                    'NGTER',     # β2 domain
                    'NLLR',      # β2 domain
                    'VTKR',      # β2 domain
                ],
                'n_terminal_signatures': [
                    'MSWKKA',    # Very specific β chain start
                    'MIHSPS',    # Alternative start
                    'SWKKA',     # Variant without Met
                ],
                'c_terminal_patterns': [
                    'VTKR',      # β2 domain
                    'RPSAE',     # C-terminal region
                    'TLLV',      # Transmembrane region
                ],
                'domain_boundaries': {
                    'beta1': (1, 90),
                    'beta2': (91, 180),
                    'beta3': (181, 210),
                },
            }
        }
    
    def _load_tcr_patterns(self) -> Dict:
        """Load TCR alpha and beta chain patterns (similar to MHC-I analyzer)."""
        return {
            'alpha': {
                'length_range': (180, 230),  # Alpha chains are typically shorter
                'optimal_length': (190, 220),
                'highly_specific_motifs': [
                    'WYQQFP',    # Very specific TCR alpha pattern
                    'WYQQ',      # Alpha-specific framework
                    'TNDYIT',    # Alpha V-region pattern
                    'QPISMDS',   # Alpha N-terminal pattern
                    'LAKTTQ',    # Specific alpha start
                ],
                'conserved_motifs': [
                    'FGXGT',    # CDR3 C-terminal
                    'WYQQKP',   # Framework regions
                    'YKFK',     # Framework 1
                    'TLTIS',    # Framework 3
                    'QLLE',     # Common N-terminal for alpha
                    'SPQFL',    # Alpha-specific pattern
                    'DSYEG',    # Alpha framework
                ],
                'n_terminal_patterns': ['QLLE', 'QVQL', 'EVQL', 'LAKTTQ'],  # Alpha N-terminal
                'c_terminal_patterns': ['FGXGT', 'FGAGT', 'FGQGT', 'FPSSP'],
            },
            'beta': {
                'length_range': (230, 290),  # Beta chains are typically longer
                'optimal_length': (240, 280),
                'highly_specific_motifs': [
                    'VKVTQS',    # Very specific beta N-terminal
                    'SSRYLV',    # Beta framework pattern
                    'LIYFSYD',   # Beta CDR2 region
                    'RTGEKV',    # Beta V-region
                    'YRQDPG',    # Beta framework
                ],
                'conserved_motifs': [
                    'FGGGT',    # CDR3 C-terminal
                    'WYQQKP',   # Framework regions
                    'MGIGV',    # Framework 1
                    'SVGD',     # Framework 2
                    'GITQ',     # Common N-terminal for beta
                    'SPKYL',    # Beta-specific pattern
                ],
                'n_terminal_patterns': ['GITQ', 'MGIT', 'GVTQ', 'VKVTQ'],  # Beta N-terminal
                'c_terminal_patterns': ['FGGGT', 'FGQGT', 'FGPGT', 'WTQDRA'],
            }
        }
    
    def _load_peptide_patterns(self) -> Dict:
        """Load peptide patterns for MHC-II (longer peptides than MHC-I)."""
        return {
            'length_range': (12, 25),  # MHC-II peptides are longer
            'optimal_length': (14, 20),  # Most common range
            'common_residues': ['A', 'L', 'V', 'I', 'F', 'Y', 'W'],  # Hydrophobic anchor residues
            'molecular_weight': (1200, 3000),  # Da
        }
    
    def analyze_pdb(self, pdb_file: str) -> Dict[str, str]:
        """
        Analyze a PDB file and identify MHC-II-TCR complex chain types.
        
        Args:
            pdb_file: Path to the PDB file
            
        Returns:
            Dictionary mapping chain IDs to chain types
        """
        if not os.path.exists(pdb_file):
            raise FileNotFoundError(f"PDB file not found: {pdb_file}")
        
        if self.verbose:
            print(f"Analyzing MHC-II complex PDB file: {pdb_file}")
        
        # Parse the structure
        structure = self.parser.get_structure('complex', pdb_file)
        
        # Extract chains and their sequences
        chain_info = self._extract_chain_info(structure)
        
        if self.verbose:
            print(f"Found {len(chain_info)} chains")
            for chain_id, info in chain_info.items():
                print(f"  Chain {chain_id}: {len(info['sequence'])} residues")
        
        # Detect multiple complexes
        complexes = self._detect_complexes(chain_info, structure)
        
        if self.verbose and len(complexes) > 1:
            print(f"\nDetected {len(complexes)} complexes:")
            for i, complex_chains in enumerate(complexes, 1):
                chain_ids = sorted(complex_chains.keys())
                print(f"  Complex {i}: Chains {', '.join(chain_ids)}")
        
        # Classify chains within each complex
        if len(complexes) > 1:
            chain_assignments = self._classify_multiple_complexes(complexes)
        else:
            # Single complex - use original logic
            chain_assignments = self._classify_chains(chain_info)
        
        if self.verbose:
            print("\nFinal chain assignments:")
            # Sort by chain ID for consistent output
            for chain_id in sorted(chain_assignments.keys()):
                print(f"  Chain {chain_id}: {chain_assignments[chain_id]}")
        
        return chain_assignments
    
    def _extract_chain_info(self, structure: Structure) -> Dict:
        """
        Extract information about each chain in the structure.
        
        Args:
            structure: BioPython structure object
            
        Returns:
            Dictionary with chain information
        """
        chain_info = {}
        
        for model in structure:
            for chain in model:
                chain_id = chain.id
                sequence = self._extract_sequence(chain)
                
                if len(sequence) > 0:  # Only include chains with sequence
                    # Calculate molecular properties
                    properties = self._calculate_properties(sequence)
                    
                    chain_info[chain_id] = {
                        'sequence': sequence,
                        'length': len(sequence),
                        'properties': properties
                    }
        
        return chain_info
    
    def _extract_sequence(self, chain) -> str:
        """Extract amino acid sequence from a chain."""
        sequence = []
        
        for residue in chain:
            if residue.has_id('CA'):  # Only standard amino acids
                res_name = residue.get_resname().strip()
                aa = self._three_to_one(res_name)
                if aa != 'X':  # Valid amino acid
                    sequence.append(aa)
        
        return ''.join(sequence)
    
    def _three_to_one(self, three_letter: str) -> str:
        """Convert three-letter amino acid code to one-letter."""
        aa_map = {
            'ALA': 'A', 'CYS': 'C', 'ASP': 'D', 'GLU': 'E', 'PHE': 'F',
            'GLY': 'G', 'HIS': 'H', 'ILE': 'I', 'LYS': 'K', 'LEU': 'L',
            'MET': 'M', 'ASN': 'N', 'PRO': 'P', 'GLN': 'Q', 'ARG': 'R',
            'SER': 'S', 'THR': 'T', 'VAL': 'V', 'TRP': 'W', 'TYR': 'Y'
        }
        return aa_map.get(three_letter, 'X')
    
    def _calculate_properties(self, sequence: str) -> Dict:
        """Calculate molecular properties of a sequence."""
        try:
            protein = ProtParam.ProteinAnalysis(sequence)
            return {
                'molecular_weight': protein.molecular_weight(),
                'gravy': protein.gravy(),  # Hydrophobicity
                'instability_index': protein.instability_index(),
                'charge': sum(1 for aa in sequence if aa in 'KR') - sum(1 for aa in sequence if aa in 'DE')
            }
        except Exception:
            return {
                'molecular_weight': len(sequence) * 110,  # Average AA weight
                'gravy': 0.0,
                'instability_index': 40.0,
                'charge': 0
            }
    
    def _detect_complexes(self, chain_info: Dict, structure) -> List[Dict]:
        """
        Detect multiple MHC-II-TCR complexes in the structure.
        Uses similar logic to MHC-I analyzer but adapted for MHC-II.
        """
        # First, try chain counting approach
        complexes_by_counting = self._detect_complexes_by_counting(chain_info)
        
        if self.verbose and len(complexes_by_counting) > 1:
            print(f"Chain counting suggests {len(complexes_by_counting)} complexes")
        
        # If counting detects multiple complexes, use that result
        if len(complexes_by_counting) > 1:
            return complexes_by_counting
        
        # If counting suggests single complex but we have many chains, try spatial clustering
        if len(chain_info) >= 8:  # Only try spatial if we have many chains
            complexes_by_spatial = self._detect_complexes_by_spatial_clustering(chain_info, structure)
            
            if self.verbose and len(complexes_by_spatial) > 1:
                print(f"Spatial clustering suggests {len(complexes_by_spatial)} complexes")
            
            return complexes_by_spatial
        
        # Default to single complex
        return complexes_by_counting
    
    def _detect_complexes_by_counting(self, chain_info: Dict) -> List[Dict]:
        """
        Detect complexes by counting potential chain types for MHC-II.
        """
        # Score all chains for all types
        potential_types = {
            'peptide': [],
            'mhc_ii_alpha': [],
            'mhc_ii_beta': [],
            'tcr_alpha': [],
            'tcr_beta': []
        }
        
        for chain_id, info in chain_info.items():
            seq = info['sequence']
            length = len(seq)
            props = info['properties']
            
            # Score for each type
            scores = {
                'peptide': self._score_peptide(seq, length, props),
                'mhc_ii_alpha': self._score_mhc_ii_alpha(seq, length),
                'mhc_ii_beta': self._score_mhc_ii_beta(seq, length),
                'tcr_alpha': self._score_tcr_alpha(seq, length),
                'tcr_beta': self._score_tcr_beta(seq, length)
            }
            
            # Add to potential types if score is reasonable
            for chain_type, score in scores.items():
                if score > 0.3:  # Threshold for being a potential candidate
                    potential_types[chain_type].append((chain_id, score))
        
        # Sort by score for each type
        for chain_type in potential_types:
            potential_types[chain_type].sort(key=lambda x: x[1], reverse=True)
        
        # Count strong candidates for unique types (no b2m in MHC-II)
        strong_peptides = len([s for _, s in potential_types['peptide'] if s > 0.6])
        strong_mhc_alpha = len([s for _, s in potential_types['mhc_ii_alpha'] if s > 0.6])
        strong_mhc_beta = len([s for _, s in potential_types['mhc_ii_beta'] if s > 0.6])
        
        # Estimate number of complexes
        num_complexes = max(strong_peptides, strong_mhc_alpha, strong_mhc_beta, 1)
        
        if num_complexes == 1:
            return [chain_info]
        
        # Group chains into complexes
        complexes = []
        assigned_chains = set()
        
        for i in range(num_complexes):
            complex_chains = {}
            
            # Assign one of each unique type per complex
            for unique_type in ['peptide', 'mhc_ii_alpha', 'mhc_ii_beta']:
                candidates = potential_types[unique_type]
                score_threshold = 0.3 if unique_type != 'mhc_ii_alpha' else 0.1  # Lower threshold for MHC
                for chain_id, score in candidates:
                    if chain_id not in assigned_chains and score > score_threshold:
                        complex_chains[chain_id] = chain_info[chain_id]
                        assigned_chains.add(chain_id)
                        break
            
            complexes.append(complex_chains)
        
        # Assign remaining TCR chains
        tcr_chains = []
        for chain_type in ['tcr_alpha', 'tcr_beta']:
            for chain_id, score in potential_types[chain_type]:
                if chain_id not in assigned_chains and score > 0.3:
                    tcr_chains.append((chain_id, chain_type, score))
        
        # Distribute TCR chains among complexes
        tcr_chains.sort(key=lambda x: x[2], reverse=True)
        
        for chain_id, chain_type, score in tcr_chains:
            # Find complex that needs this type
            best_complex_idx = 0
            for i, complex_chains in enumerate(complexes):
                # Check how many TCR chains this complex already has
                tcr_count = sum(1 for cid in complex_chains.keys() 
                              if len(chain_info[cid]['sequence']) > 150)
                if tcr_count < 2:  # Each complex should have 2 TCR chains
                    best_complex_idx = i
                    break
            
            complexes[best_complex_idx][chain_id] = chain_info[chain_id]
            assigned_chains.add(chain_id)
        
        # Add any remaining unassigned chains to complexes
        for chain_id in chain_info:
            if chain_id not in assigned_chains:
                # Add to the smallest complex
                smallest_idx = min(range(len(complexes)), key=lambda i: len(complexes[i]))
                complexes[smallest_idx][chain_id] = chain_info[chain_id]
        
        # Filter out empty complexes
        complexes = [c for c in complexes if c]
        
        return complexes if len(complexes) > 1 else [chain_info]
    
    def _detect_complexes_by_spatial_clustering(self, chain_info: Dict, structure) -> List[Dict]:
        """
        Detect complexes by spatial clustering of chain centers.
        (Similar to MHC-I analyzer)
        """
        try:
            # Calculate center of mass for each chain
            chain_centers = {}
            
            for model in structure:
                for chain in model:
                    chain_id = chain.id
                    if chain_id in chain_info:
                        # Calculate center of mass
                        coords = []
                        for residue in chain:
                            if residue.has_id('CA'):  # Use alpha carbon
                                ca_atom = residue['CA']
                                coords.append(ca_atom.coord)
                        
                        if coords:
                            center = np.mean(coords, axis=0)
                            chain_centers[chain_id] = center
            
            if len(chain_centers) < 2:
                return [chain_info]
            
            # Calculate pairwise distances
            chain_ids = list(chain_centers.keys())
            distances = {}
            
            for i, chain1 in enumerate(chain_ids):
                for j, chain2 in enumerate(chain_ids[i+1:], i+1):
                    dist = np.linalg.norm(chain_centers[chain1] - chain_centers[chain2])
                    distances[(chain1, chain2)] = dist
            
            # Simple clustering: chains within 50Å are in the same complex
            cluster_threshold = 50.0
            clusters = []
            assigned = set()
            
            for chain_id in chain_ids:
                if chain_id in assigned:
                    continue
                
                cluster = {chain_id}
                to_check = [chain_id]
                
                while to_check:
                    current = to_check.pop()
                    for other_chain in chain_ids:
                        if other_chain in cluster or other_chain in assigned:
                            continue
                        
                        # Check distance
                        pair = tuple(sorted([current, other_chain]))
                        if pair in distances and distances[pair] < cluster_threshold:
                            cluster.add(other_chain)
                            to_check.append(other_chain)
                
                clusters.append(cluster)
                assigned.update(cluster)
            
            # Convert clusters to complex format
            complexes = []
            for cluster in clusters:
                complex_chains = {cid: chain_info[cid] for cid in cluster}
                complexes.append(complex_chains)
            
            return complexes if len(complexes) > 1 else [chain_info]
        
        except Exception as e:
            if self.verbose:
                print(f"Spatial clustering failed: {e}, falling back to single complex")
            return [chain_info]
    
    def _classify_chains(self, chain_info: Dict) -> Dict[str, str]:
        """
        Classify chains in a single MHC-II complex.
        """
        assignments = {}
        
        # Two-pass classification system
        # Pass 1: Calculate all scores
        all_scores = {}
        for chain_id, info in chain_info.items():
            seq = info['sequence']
            length = len(seq)
            props = info['properties']
            
            scores = {
                'peptide': self._score_peptide(seq, length, props),
                'mhc_ii_alpha': self._score_mhc_ii_alpha(seq, length),
                'mhc_ii_beta': self._score_mhc_ii_beta(seq, length),
                'tcr_alpha': self._score_tcr_alpha(seq, length),
                'tcr_beta': self._score_tcr_beta(seq, length)
            }
            all_scores[chain_id] = scores
        
        # Pass 2: Assign chains based on best scores, avoiding conflicts
        used_unique_types = set()
        
        # Sort chains by length (shortest first for peptide priority)
        sorted_chains = sorted(chain_info.items(), key=lambda x: x[1]['length'])
        
        for chain_id, info in sorted_chains:
            scores = all_scores[chain_id]
            
            # Find best scoring type that hasn't been used (for unique types)
            best_type = max(scores.items(), key=lambda x: x[1])
            best_score = best_type[1]
            best_chain_type = best_type[0]
            
            # Check if this is a unique type already assigned
            unique_types = {'peptide', 'mhc_ii_alpha', 'mhc_ii_beta'}
            
            if best_chain_type in unique_types and best_chain_type in used_unique_types:
                # Find next best type that's not unique or not used
                sorted_scores = sorted(scores.items(), key=lambda x: x[1], reverse=True)
                assigned = False
                for chain_type, score in sorted_scores:
                    if score > 0.3:  # Minimum threshold
                        if chain_type not in unique_types or chain_type not in used_unique_types:
                            assignments[chain_id] = chain_type
                            if chain_type in unique_types:
                                used_unique_types.add(chain_type)
                            assigned = True
                            break
                
                if not assigned:
                    assignments[chain_id] = 'unknown'
            else:
                # Can assign the best type
                if best_score > 0.3:
                    assignments[chain_id] = best_chain_type
                    if best_chain_type in unique_types:
                        used_unique_types.add(best_chain_type)
                else:
                    assignments[chain_id] = 'unknown'
        
        return assignments
    
    def _classify_multiple_complexes(self, complexes: List[Dict]) -> Dict[str, str]:
        """
        Classify chains within multiple MHC-II complexes.
        """
        all_assignments = {}
        
        for i, complex_chains in enumerate(complexes):
            if self.verbose:
                print(f"\nClassifying MHC-II Complex {i+1} ({len(complex_chains)} chains):")
                for chain_id, info in complex_chains.items():
                    print(f"  Chain {chain_id}: {len(info['sequence'])} residues")
            
            # Classify this complex
            complex_assignments = self._classify_chains(complex_chains)
            
            # Add complex number to chain assignments only if multiple complexes
            for chain_id, chain_type in complex_assignments.items():
                if len(complexes) > 1:
                    # Multiple complexes - add complex number
                    if chain_type != 'unknown':
                        all_assignments[chain_id] = f"{chain_type}_complex{i+1}"
                    else:
                        all_assignments[chain_id] = 'unknown'
                else:
                    # Single complex - use original names
                    all_assignments[chain_id] = chain_type
        
        return all_assignments
    
    # Scoring functions for each chain type
    
    def _score_peptide(self, sequence: str, length: int, properties: Dict) -> float:
        """Score likelihood of being a peptide antigen (MHC-II version - longer)."""
        score = 0.0
        
        # Length scoring for MHC-II peptides (longer than MHC-I)
        peptide_range = self.peptide_patterns['length_range']
        optimal_range = self.peptide_patterns['optimal_length']
        
        if peptide_range[0] <= length <= peptide_range[1]:
            score += 0.8  # Strong preference for correct length
            if optimal_range[0] <= length <= optimal_range[1]:
                score += 0.2  # Bonus for optimal length
        elif length < peptide_range[0]:
            # Too short for MHC-II peptide
            score -= 0.5
        elif length > peptide_range[1] + 5:
            # Too long for peptide
            score -= 0.8
        
        # Molecular weight check
        mw_range = self.peptide_patterns['molecular_weight']
        mw = properties.get('molecular_weight', 0)
        if mw_range[0] <= mw <= mw_range[1]:
            score += 0.2
        
        # Hydrophobic residues (common in peptide anchors)
        hydrophobic_count = sum(1 for aa in sequence if aa in self.peptide_patterns['common_residues'])
        hydrophobic_fraction = hydrophobic_count / length if length > 0 else 0
        score += hydrophobic_fraction * 0.3
        
        # Very short sequences are likely peptides
        if length <= 11:
            score += 0.4
        
        return min(max(score, 0.0), 1.0)
    
    def _score_mhc_ii_alpha(self, sequence: str, length: int) -> float:
        """Score likelihood of being MHC-II alpha chain."""
        score = 0.0
        
        # Length scoring
        alpha_range = self.mhc_ii_patterns['alpha']['length_range']
        optimal_range = self.mhc_ii_patterns['alpha']['optimal_length']
        
        if alpha_range[0] <= length <= alpha_range[1]:
            score += 0.7  # Strong weight for correct length
            if optimal_range[0] <= length <= optimal_range[1]:
                score += 0.3  # Bonus for optimal length
        elif length < alpha_range[0] - 10 or length > alpha_range[1] + 20:
            score -= 0.6  # Strong penalty for wrong length
        
        # Check for highly specific MHC-II α patterns
        n_terminal = sequence[:50] if len(sequence) >= 50 else sequence
        c_terminal = sequence[-50:] if len(sequence) >= 50 else sequence
        
        # Highly specific alpha motifs
        alpha_specific_found = 0
        for pattern in self.mhc_ii_patterns['alpha']['highly_specific_motifs']:
            if pattern in sequence:
                alpha_specific_found += 1
                score += 0.4  # High score for specific patterns
        
        # Multiple specific patterns = strong evidence
        if alpha_specific_found >= 2:
            score += 0.5
        
        # N-terminal signature patterns
        for pattern in self.mhc_ii_patterns['alpha']['n_terminal_signatures']:
            if pattern in n_terminal:
                score += 0.5  # Very strong MHC-II α indicator
                break
        
        # C-terminal patterns
        for pattern in self.mhc_ii_patterns['alpha']['c_terminal_patterns']:
            if pattern in c_terminal:
                score += 0.2
                break
        
        # Regular conserved motifs
        motifs = self.mhc_ii_patterns['alpha']['conserved_motifs']
        motif_count = sum(1 for motif in motifs if motif in sequence)
        score += motif_count * 0.1
        
        # Penalty for TCR-like patterns
        tcr_patterns = (self.tcr_patterns['alpha']['highly_specific_motifs'] + 
                       self.tcr_patterns['beta']['highly_specific_motifs'])
        if any(pattern in sequence for pattern in tcr_patterns):
            score -= 0.6
        
        return min(max(score, 0.0), 1.0)
    
    def _score_mhc_ii_beta(self, sequence: str, length: int) -> float:
        """Score likelihood of being MHC-II beta chain."""
        score = 0.0
        
        # Length scoring
        beta_range = self.mhc_ii_patterns['beta']['length_range']
        optimal_range = self.mhc_ii_patterns['beta']['optimal_length']
        
        if beta_range[0] <= length <= beta_range[1]:
            score += 0.7  # Strong weight for correct length
            if optimal_range[0] <= length <= optimal_range[1]:
                score += 0.3  # Bonus for optimal length
        elif length < beta_range[0] - 10 or length > beta_range[1] + 20:
            score -= 0.6  # Strong penalty for wrong length
        
        # Check for highly specific MHC-II β patterns
        n_terminal = sequence[:50] if len(sequence) >= 50 else sequence
        c_terminal = sequence[-50:] if len(sequence) >= 50 else sequence
        
        # Highly specific beta motifs
        beta_specific_found = 0
        for pattern in self.mhc_ii_patterns['beta']['highly_specific_motifs']:
            if pattern in sequence:
                beta_specific_found += 1
                score += 0.4  # High score for specific patterns
        
        # Multiple specific patterns = strong evidence
        if beta_specific_found >= 2:
            score += 0.5
        
        # N-terminal signature patterns
        for pattern in self.mhc_ii_patterns['beta']['n_terminal_signatures']:
            if pattern in n_terminal:
                score += 0.5  # Very strong MHC-II β indicator
                break
        
        # C-terminal patterns
        for pattern in self.mhc_ii_patterns['beta']['c_terminal_patterns']:
            if pattern in c_terminal:
                score += 0.2
                break
        
        # Regular conserved motifs
        motifs = self.mhc_ii_patterns['beta']['conserved_motifs']
        motif_count = sum(1 for motif in motifs if motif in sequence)
        score += motif_count * 0.1
        
        # Penalty for TCR-like patterns
        tcr_patterns = (self.tcr_patterns['alpha']['highly_specific_motifs'] + 
                       self.tcr_patterns['beta']['highly_specific_motifs'])
        if any(pattern in sequence for pattern in tcr_patterns):
            score -= 0.6
        
        return min(max(score, 0.0), 1.0)
    
    def _score_tcr_alpha(self, sequence: str, length: int) -> float:
        """Score likelihood of being TCR alpha chain (same as MHC-I analyzer)."""
        score = 0.0
        
        # Length check with stronger discrimination
        alpha_range = self.tcr_patterns['alpha']['length_range']
        optimal_range = self.tcr_patterns['alpha']['optimal_length']
        
        if alpha_range[0] <= length <= alpha_range[1]:
            score += 0.5  # Strong length preference
            # Bonus for optimal length range
            if optimal_range[0] <= length <= optimal_range[1]:
                score += 0.3
        elif length > alpha_range[1]:
            # Penalty for being too long (likely beta)
            score -= 0.4
        elif length < alpha_range[0]:
            # Penalty for being too short
            score -= 0.3
        
        # Check for highly specific alpha patterns
        n_terminal = sequence[:30] if len(sequence) >= 30 else sequence
        c_terminal = sequence[-30:] if len(sequence) >= 30 else sequence
        
        # Highly specific alpha motifs
        alpha_specific_found = 0
        for pattern in self.tcr_patterns['alpha']['highly_specific_motifs']:
            if pattern in sequence:
                alpha_specific_found += 1
                score += 0.4  # High score for specific patterns
        
        # Multiple specific patterns = strong evidence
        if alpha_specific_found >= 2:
            score += 0.4
        
        # N-terminal patterns (very distinctive for alpha)
        n_patterns = self.tcr_patterns['alpha']['n_terminal_patterns']
        for pattern in n_patterns:
            if pattern in n_terminal:
                score += 0.35
                break
        
        # C-terminal patterns
        c_patterns = self.tcr_patterns['alpha']['c_terminal_patterns']
        for pattern in c_patterns:
            if pattern in c_terminal:
                score += 0.25
                break
        
        # Regular motif checks
        motifs = self.tcr_patterns['alpha']['conserved_motifs']
        motif_count = sum(1 for motif in motifs if motif in sequence)
        score += motif_count * 0.1
        
        # Strong preference for shorter chains (alpha characteristic)
        if length < 220:
            score += 0.2
        elif length > 250:
            score -= 0.5  # Strong penalty for long chains
        
        # Penalty for MHC-II-specific patterns
        mhc_alpha_patterns = self.mhc_ii_patterns['alpha']['highly_specific_motifs']
        mhc_beta_patterns = self.mhc_ii_patterns['beta']['highly_specific_motifs']
        if any(pattern in sequence for pattern in mhc_alpha_patterns + mhc_beta_patterns):
            score -= 0.8  # Strong penalty for MHC patterns
        
        # Penalty for beta-specific patterns
        beta_specific = self.tcr_patterns['beta']['highly_specific_motifs']
        if any(pattern in sequence for pattern in beta_specific):
            score -= 0.6  # Penalty for beta-specific patterns
        
        return min(max(score, 0.0), 1.0)
    
    def _score_tcr_beta(self, sequence: str, length: int) -> float:
        """Score likelihood of being TCR beta chain (same as MHC-I analyzer)."""
        score = 0.0
        
        # Length check with stronger discrimination
        beta_range = self.tcr_patterns['beta']['length_range']
        optimal_range = self.tcr_patterns['beta']['optimal_length']
        
        if beta_range[0] <= length <= beta_range[1]:
            score += 0.5  # Strong length preference
            # Bonus for optimal length range
            if optimal_range[0] <= length <= optimal_range[1]:
                score += 0.3
        elif length < beta_range[0]:
            # Penalty for being too short (likely alpha)
            score -= 0.4
        elif length > beta_range[1]:
            # Penalty for being too long (likely MHC)
            score -= 0.3
        
        # Check for highly specific beta patterns
        n_terminal = sequence[:30] if len(sequence) >= 30 else sequence
        c_terminal = sequence[-30:] if len(sequence) >= 30 else sequence
        
        # Highly specific beta motifs
        beta_specific_found = 0
        for pattern in self.tcr_patterns['beta']['highly_specific_motifs']:
            if pattern in sequence:
                beta_specific_found += 1
                score += 0.4  # High score for specific patterns
        
        # Multiple specific patterns = strong evidence
        if beta_specific_found >= 2:
            score += 0.4
        
        # N-terminal patterns (very distinctive for beta)
        n_patterns = self.tcr_patterns['beta']['n_terminal_patterns']
        for pattern in n_patterns:
            if pattern in n_terminal:
                score += 0.35
                break
        
        # C-terminal patterns
        c_patterns = self.tcr_patterns['beta']['c_terminal_patterns']
        for pattern in c_patterns:
            if pattern in c_terminal:
                score += 0.25
                break
        
        # Regular motif checks
        motifs = self.tcr_patterns['beta']['conserved_motifs']
        motif_count = sum(1 for motif in motifs if motif in sequence)
        score += motif_count * 0.1
        
        # Strong preference for longer chains (beta characteristic)
        if length > 230 and length < 290:
            score += 0.2
        elif length < 220:
            score -= 0.5  # Strong penalty for short chains
        elif length > 300:
            score -= 0.6  # Very strong penalty for very long chains (likely MHC)
        
        # Penalty for MHC-II-specific patterns
        mhc_alpha_patterns = self.mhc_ii_patterns['alpha']['highly_specific_motifs']
        mhc_beta_patterns = self.mhc_ii_patterns['beta']['highly_specific_motifs']
        if any(pattern in sequence for pattern in mhc_alpha_patterns + mhc_beta_patterns):
            score -= 0.8  # Strong penalty for MHC patterns
        
        # Penalty for alpha-specific patterns
        alpha_specific = self.tcr_patterns['alpha']['highly_specific_motifs']
        if any(pattern in sequence for pattern in alpha_specific):
            score -= 0.6  # Penalty for alpha-specific patterns
        
        return min(max(score, 0.0), 1.0)
    
    def save_analysis_report(self, assignments: Dict[str, str], 
                           output_file: str, chain_info: Dict = None):
        """
        Save MHC-II analysis report to file.
        
        Args:
            assignments: Chain assignments
            output_file: Output file path
            chain_info: Additional chain information
        """
        with open(output_file, 'w') as f:
            f.write("MHC-II-TCR Complex Analysis Report\n")
            f.write("=" * 40 + "\n\n")
            
            # Sort by chain ID for consistent output
            for chain_id in sorted(assignments.keys()):
                chain_type = assignments[chain_id]
                f.write(f"Chain {chain_id}: {chain_type}\n")
                
                if chain_info and chain_id in chain_info:
                    info = chain_info[chain_id]
                    f.write(f"  Length: {info['length']} residues\n")
                    f.write(f"  Molecular Weight: {info['properties'].get('molecular_weight', 'N/A'):.1f} Da\n")
                    f.write(f"  Sequence: {info['sequence'][:50]}...\n\n")
        
        if self.verbose:
            print(f"Analysis report saved to: {output_file}")


def main():
    """Command-line interface for MHC-II analyzer."""
    import sys
    
    if len(sys.argv) != 2:
        print("Usage: python mhc_ii_analyzer.py <pdb_file>")
        sys.exit(1)
    
    pdb_file = sys.argv[1]
    
    try:
        analyzer = pMHCIITCRAnalyzer(verbose=True)
        assignments = analyzer.analyze_pdb(pdb_file)
        
        print("\nFinal Chain Assignments:")
        print("-" * 30)
        for chain_id in sorted(assignments.keys()):
            print(f"Chain {chain_id}: {assignments[chain_id]}")
        
    except Exception as e:
        print(f"Error: {e}")
        sys.exit(1)


if __name__ == "__main__":
    main()