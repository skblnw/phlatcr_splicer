#!/usr/bin/env python3
"""
pHLA-TCR Complex Structure Analyzer

A Python tool for analyzing protein complex structure PDB files, specifically 
designed for pHLA-TCR complexes. This tool identifies and classifies protein 
chains even when chain IDs are inconsistent.

Author: [Your Name]
License: BSD 3-Clause
"""

import os
import sys
import re
import warnings
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
    sys.exit(1)


class pMHCITCRAnalyzer:
    """
    Main class for analyzing pHLA-TCR complex structures.
    """
    
    def __init__(self, verbose: bool = True):
        """
        Initialize the analyzer.
        
        Args:
            verbose: Whether to print verbose output
        """
        self.verbose = verbose
        self.parser = PDBParser(QUIET=not verbose)
        
        # Define known sequence patterns and characteristics
        self.mhc_patterns = self._load_mhc_patterns()
        self.tcr_patterns = self._load_tcr_patterns()
        self.b2m_patterns = self._load_b2m_patterns()
        
    def _load_mhc_patterns(self) -> Dict:
        """Load MHC class I heavy chain patterns and characteristics."""
        return {
            'length_range': (230, 380),  # Extended range for shorter variants
            'optimal_length': (265, 285),  # Most common MHC heavy chain lengths
            'acceptable_length': (230, 300),  # Acceptable range for most MHC
            'highly_specific_motifs': [
                'GSHSMR',    # Very specific MHC Heavy N-terminal
                'GSHSM',     # Highly specific MHC marker
                'YFYTS',     # Alpha1 domain signature
                'YFYT',      # Variant of above
                'MRYFYT',    # Extended N-terminal pattern
                'RYFYTS',    # Extended pattern
            ],
            'conserved_motifs': [
                'WSDRVI',   # Alpha2 domain motif
                'FKAFLKQ',  # Alpha3 domain motif
                'TGAASC',   # Common motif
                'RGEC',     # Additional pattern
                'ASSPRE',   # Found in many MHC sequences
                'EPRAP',    # N-terminal region pattern
            ],
            'n_terminal_signatures': [
                'GSHSMR',    # Most specific
                'GSHSM',     # Highly specific
                'SHSMR',     # Variant
                'HSMRY',     # Variant
            ],
            'c_terminal_patterns': [
                'KWAAV',     # C-terminal pattern
                'WAAVV',     # C-terminal variant
                'TQIVSAE',   # C-terminal signature
            ],
            'anti_tcr_patterns': [  # Patterns that rule out TCR
                'GSHSM',     # Never in TCR
                'YFYTS',     # Never in TCR
                'EPRAP',     # MHC-specific
            ],
            'domain_boundaries': {
                'alpha1': (1, 90),
                'alpha2': (91, 180),
                'alpha3': (181, 276),
            },
            'disulfide_cysteines': [25, 101, 164, 236],  # Typical cysteine positions
            'hydrophobic_residues_high': True,  # MHC heavy chains are typically hydrophobic
        }
    
    def _load_tcr_patterns(self) -> Dict:
        """Load TCR alpha and beta chain patterns."""
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
                'v_gene_patterns': ['TRAV', 'TCRAV'],
                'j_gene_patterns': ['TRAJ', 'TCRAJ'],
                'cdr3_patterns': ['CAS', 'CAV', 'CAI'],  # Common CDR3 starts
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
                'v_gene_patterns': ['TRBV', 'TCRBV'],
                'j_gene_patterns': ['TRBJ', 'TCRBJ'],
                'cdr3_patterns': ['CAS', 'CAW', 'CAT'],  # Common CDR3 starts
                'n_terminal_patterns': ['GITQ', 'MGIT', 'GVTQ', 'VKVTQ'],  # Beta N-terminal
                'c_terminal_patterns': ['FGGGT', 'FGQGT', 'FGPGT', 'WTQDRA'],
            }
        }
    
    def _load_b2m_patterns(self) -> Dict:
        """Load β2-microglobulin patterns."""
        return {
            'length_range': (95, 105),  # Very consistent length
            'conserved_motifs': [
                'IQRTPK',   # N-terminal
                'VNHVTL',   # Middle region
                'SRNLTKDR', # C-terminal
            ],
            'molecular_weight': (11000, 13000),  # Da
        }
    
    def analyze_pdb(self, pdb_file: str) -> Dict[str, str]:
        """
        Analyze a PDB file and identify chain types.
        
        Args:
            pdb_file: Path to the PDB file
            
        Returns:
            Dictionary mapping chain IDs to chain types
        """
        if not os.path.exists(pdb_file):
            raise FileNotFoundError(f"PDB file not found: {pdb_file}")
        
        if self.verbose:
            print(f"Analyzing PDB file: {pdb_file}")
        
        # Parse the structure
        structure = self.parser.get_structure('complex', pdb_file)
        
        # Extract chains and their sequences
        chain_info = self._extract_chain_info(structure)
        
        if self.verbose:
            print(f"Found {len(chain_info)} chains")
            for chain_id, info in chain_info.items():
                print(f"  Chain {chain_id}: {len(info['sequence'])} residues")
        
        # Try to detect complexes from REMARK hints first (e.g., REMARK 350 chains groups)
        complexes = self._detect_complexes_from_remarks(pdb_file, chain_info)
        
        # If REMARKs did not provide grouping, fall back to existing heuristics
        if len(complexes) <= 1:
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

    def _detect_complexes_from_remarks(self, pdb_file: str, chain_info: Dict) -> List[Dict]:
        """
        Use REMARK lines (especially REMARK 350) to pre-group chains into complexes.
        
        Many multi-complex pHLA-TCR PDBs list biological assembly groups as
        "APPLY THE FOLLOWING TO CHAINS: A, B, C, M, P". These groups typically
        correspond to one pHLA-TCR complex (5 chains: MHC heavy, β2m, peptide,
        TCR α, TCR β).
        
        Args:
            pdb_file: Path to the PDB file
            chain_info: Extracted chain information
        
        Returns:
            List of dictionaries, each containing the chains for one complex.
            Returns [chain_info] if no usable REMARK grouping is found.
        """
        try:
            groups: List[List[str]] = []
            with open(pdb_file, 'r') as fh:
                for line in fh:
                    # Focus on REMARK 350 lines that specify chain application
                    if line.startswith('REMARK 350') and 'APPLY THE FOLLOWING TO CHAINS:' in line:
                        # Extract substring after the colon
                        parts = line.split('APPLY THE FOLLOWING TO CHAINS:')
                        if len(parts) < 2:
                            continue
                        tail = parts[1].strip()
                        # Remove trailing whitespace and anything after two or more spaces (often BIOMT follows on later lines)
                        tail = tail.split('  ')[0].strip()
                        # Split by commas and spaces to get chain ids
                        raw_tokens = [t.strip() for t in tail.split(',') if t.strip()]
                        # Keep only first character of token (classic PDB chain IDs are single char)
                        chains = [tok[0] for tok in raw_tokens if len(tok) > 0]
                        # Basic sanity: groups that look like plausible complexes (4-6 chains)
                        if len(chains) >= 4:
                            groups.append(chains)
            # Deduplicate identical groups
            unique_groups: List[List[str]] = []
            seen = set()
            for grp in groups:
                key = tuple(sorted(grp))
                if key not in seen:
                    seen.add(key)
                    unique_groups.append(grp)
            # Map groups to chain_info
            if len(unique_groups) >= 2:
                complexes: List[Dict] = []
                for grp in unique_groups:
                    complex_chains: Dict = {}
                    for cid in grp:
                        if cid in chain_info:
                            complex_chains[cid] = chain_info[cid]
                    # Only include non-empty groups
                    if complex_chains:
                        complexes.append(complex_chains)
                # Ensure we truly have multiple non-empty complexes
                if len([c for c in complexes if len(c) >= 2]) >= 2:
                    return complexes
            return [chain_info]
        except Exception:
            return [chain_info]
    
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
                chain_id = chain.get_id()
                
                # Extract sequence
                sequence = self._get_chain_sequence(chain)
                
                # Calculate basic properties
                properties = self._calculate_chain_properties(sequence, chain)
                
                chain_info[chain_id] = {
                    'sequence': sequence,
                    'chain_object': chain,
                    'length': len(sequence),
                    'properties': properties
                }
        
        return chain_info
    
    def _get_chain_sequence(self, chain: Chain) -> str:
        """
        Extract amino acid sequence from a chain.
        
        Args:
            chain: BioPython chain object
            
        Returns:
            Amino acid sequence as string
        """
        sequence = ""
        
        for residue in chain:
            if residue.get_id()[0] == ' ':  # Standard amino acid
                resname = residue.get_resname()
                # Convert 3-letter to 1-letter code
                aa = self._three_to_one(resname)
                if aa:
                    sequence += aa
        
        return sequence
    
    def _three_to_one(self, three_letter: str) -> str:
        """Convert 3-letter amino acid code to 1-letter."""
        conversion = {
            'ALA': 'A', 'ARG': 'R', 'ASN': 'N', 'ASP': 'D', 'CYS': 'C',
            'GLU': 'E', 'GLN': 'Q', 'GLY': 'G', 'HIS': 'H', 'ILE': 'I',
            'LEU': 'L', 'LYS': 'K', 'MET': 'M', 'PHE': 'F', 'PRO': 'P',
            'SER': 'S', 'THR': 'T', 'TRP': 'W', 'TYR': 'Y', 'VAL': 'V'
        }
        return conversion.get(three_letter.upper(), '')
    
    def _calculate_chain_properties(self, sequence: str, chain: Chain) -> Dict:
        """
        Calculate various properties of a protein chain.
        
        Args:
            sequence: Amino acid sequence
            chain: BioPython chain object
            
        Returns:
            Dictionary of calculated properties
        """
        properties = {}
        
        if sequence:
            # Calculate molecular weight and other properties
            try:
                protein_analysis = ProtParam.ProteinAnalysis(sequence)
                properties['molecular_weight'] = protein_analysis.molecular_weight()
                properties['aromaticity'] = protein_analysis.aromaticity()
                properties['instability_index'] = protein_analysis.instability_index()
                properties['isoelectric_point'] = protein_analysis.isoelectric_point()
            except:
                properties['molecular_weight'] = 0
                properties['aromaticity'] = 0
                properties['instability_index'] = 0
                properties['isoelectric_point'] = 7.0
        
        # Count secondary structure elements (simplified)
        properties['num_residues'] = len(list(chain.get_residues()))
        
        return properties
    
    def _classify_chains(self, chain_info: Dict) -> Dict[str, str]:
        """
        Classify chains based on sequence and structural features.
        
        Args:
            chain_info: Dictionary with chain information
            
        Returns:
            Dictionary mapping chain IDs to chain types
        """
        assignments = {}
        used_types = set()
        
        # Sort chains by length for better classification
        sorted_chains = sorted(chain_info.items(), 
                             key=lambda x: x[1]['length'], reverse=True)
        
        # First pass: calculate scores for all chain types
        chain_scores = {}
        for chain_id, info in sorted_chains:
            sequence = info['sequence']
            length = info['length']
            properties = info['properties']
            
            scores = {}
            scores['peptide'] = self._score_peptide(sequence, length, properties)
            scores['b2m'] = self._score_b2m(sequence, length, properties)
            scores['mhc_heavy'] = self._score_mhc(sequence, length, properties)
            scores['tcr_alpha'] = self._score_tcr_alpha(sequence, length)
            scores['tcr_beta'] = self._score_tcr_beta(sequence, length)
            
            chain_scores[chain_id] = scores
        
        # Second pass: assign based on best unique matches
        used_types = set()
        for chain_id, info in sorted_chains:
            sequence = info['sequence']
            length = info['length']
            properties = info['properties']
            
            # Get scores for this chain
            scores = chain_scores[chain_id]
            
            # Find the best available type
            available_scores = {ctype: score for ctype, score in scores.items() 
                              if ctype not in used_types}
            
            if available_scores:
                best_type = max(available_scores.items(), key=lambda x: x[1])
                if best_type[1] > 0.3:  # Minimum confidence threshold
                    assignments[chain_id] = best_type[0]
                    used_types.add(best_type[0])
                else:
                    # Fallback to length-based assignment
                    chain_type = self._fallback_assignment(length, used_types)
                    assignments[chain_id] = chain_type
                    used_types.add(chain_type)
            else:
                assignments[chain_id] = 'unknown'
        
        # Post-process to ensure we have reasonable assignments
        assignments = self._validate_assignments(assignments, chain_info)
        
        return assignments
    
    def _determine_chain_type(self, sequence: str, length: int, 
                            properties: Dict, used_types: Set) -> str:
        """
        Determine the type of a single chain using comprehensive scoring.
        
        Args:
            sequence: Amino acid sequence
            length: Sequence length
            properties: Chain properties
            used_types: Already assigned chain types
            
        Returns:
            Predicted chain type
        """
        # Calculate scores for each possible chain type
        scores = {}
        
        if 'peptide' not in used_types:
            scores['peptide'] = self._score_peptide(sequence, length, properties)
        
        if 'b2m' not in used_types:
            scores['b2m'] = self._score_b2m(sequence, length, properties)
        
        if 'mhc_heavy' not in used_types:
            scores['mhc_heavy'] = self._score_mhc(sequence, length, properties)
        
        if 'tcr_alpha' not in used_types:
            scores['tcr_alpha'] = self._score_tcr_alpha(sequence, length)
        
        if 'tcr_beta' not in used_types:
            scores['tcr_beta'] = self._score_tcr_beta(sequence, length)
        
        # Find the highest scoring type
        if scores:
            best_type = max(scores.items(), key=lambda x: x[1])
            if best_type[1] > 0.3:  # Minimum confidence threshold
                return best_type[0]
        
        # Fallback to length-based assignment
        return self._fallback_assignment(length, used_types)
    
    def _score_peptide(self, sequence: str, length: int, properties: Dict) -> float:
        """Score likelihood of being a peptide antigen."""
        score = 0.0
        
        # Length scoring - peptides are typically 8-15 residues
        if 8 <= length <= 15:
            score += 0.8
        elif length < 8:
            score += 0.4
        elif length <= 20:
            score += 0.2
        else:
            score -= 0.3
        
        # Very short sequences are likely peptides
        if length < 12:
            score += 0.3
        
        return min(max(score, 0.0), 1.0)
    
    def _score_b2m(self, sequence: str, length: int, properties: Dict) -> float:
        """Score likelihood of being β2-microglobulin."""
        score = 0.0
        
        # Length check - β2m is very consistent around 99 residues
        b2m_range = self.b2m_patterns['length_range']
        if b2m_range[0] <= length <= b2m_range[1]:
            score += 0.5
            # Bonus for being close to 99
            if abs(length - 99) <= 3:
                score += 0.2
        
        # Check conserved motifs
        motif_score = sum(1 for motif in self.b2m_patterns['conserved_motifs'] 
                         if motif in sequence)
        score += motif_score * 0.15
        
        # Molecular weight check
        mw = properties.get('molecular_weight', 0)
        if 11000 <= mw <= 13000:
            score += 0.15
        
        return min(score, 1.0)
    
    def _score_mhc(self, sequence: str, length: int, properties: Dict) -> float:
        """Score likelihood of being MHC heavy chain."""
        score = 0.0
        
        # Strong length preference for MHC
        mhc_range = self.mhc_patterns['length_range']
        optimal_range = self.mhc_patterns['optimal_length']
        
        if mhc_range[0] <= length <= mhc_range[1]:
            score += 0.7  # Strong weight for correct length
            # Extra bonus for optimal MHC length
            if optimal_range[0] <= length <= optimal_range[1]:
                score += 0.4  # Very strong bonus
            # Acceptable bonus for reasonable MHC length
            elif self.mhc_patterns['acceptable_length'][0] <= length <= self.mhc_patterns['acceptable_length'][1]:
                score += 0.2
        elif length < 200:
            # Strong penalty for being too short (can't be MHC)
            score -= 0.8
        
        # HIGHLY SPECIFIC MHC PATTERNS - These are almost diagnostic
        n_terminal = sequence[:50] if len(sequence) >= 50 else sequence
        c_terminal = sequence[-50:] if len(sequence) >= 50 else sequence
        
        # Check for highly specific MHC signatures
        highly_specific_found = 0
        for pattern in self.mhc_patterns['highly_specific_motifs']:
            if pattern in sequence:
                highly_specific_found += 1
                score += 0.5  # Very high score for each specific pattern
        
        # If we find 2+ highly specific patterns, this is almost certainly MHC
        if highly_specific_found >= 2:
            score += 0.6  # Massive bonus for multiple specific patterns
        # Even one highly specific pattern is very strong evidence
        elif highly_specific_found >= 1:
            score += 0.4  # Strong bonus for at least one specific pattern
        
        # N-terminal signature patterns (very diagnostic)
        n_signature_found = False
        for pattern in self.mhc_patterns['n_terminal_signatures']:
            if pattern in n_terminal:
                score += 0.6  # Very strong MHC indicator
                n_signature_found = True
                break
        
        # Check conserved motifs
        motif_score = sum(1 for motif in self.mhc_patterns['conserved_motifs'] 
                         if motif in sequence)
        score += motif_score * 0.1
        
        # MHC heavy chains should NOT have TCR patterns
        # Penalty for TCR-like patterns
        tcr_beta_patterns = self.tcr_patterns['beta']['n_terminal_patterns']
        tcr_alpha_patterns = self.tcr_patterns['alpha']['n_terminal_patterns']
        
        if any(pattern in n_terminal for pattern in tcr_beta_patterns + tcr_alpha_patterns):
            score -= 0.4  # Strong penalty for TCR patterns
        
        # Check for cysteine residues at expected positions
        cys_score = 0
        for pos in self.mhc_patterns['disulfide_cysteines']:
            if pos < len(sequence) and sequence[pos-1] == 'C':  # 1-indexed to 0-indexed
                cys_score += 1
        score += (cys_score / len(self.mhc_patterns['disulfide_cysteines'])) * 0.2
        
        # MHC chains are typically the longest protein in pHLA-TCR complexes
        if length > 270:
            score += 0.2
        
        return min(max(score, 0.0), 1.0)
    
    def _fallback_assignment(self, length: int, used_types: Set) -> str:
        """Fallback assignment based on length when scoring fails."""
        if 'peptide' not in used_types and length < 20:
            return 'peptide'
        elif 'b2m' not in used_types and length < 150:
            return 'b2m'
        elif 'mhc_heavy' not in used_types and length > 250:
            return 'mhc_heavy'
        elif 'tcr_alpha' not in used_types:
            return 'tcr_alpha'
        elif 'tcr_beta' not in used_types:
            return 'tcr_beta'
        else:
            return 'unknown'
    
    def _score_tcr_alpha(self, sequence: str, length: int) -> float:
        """Score likelihood of being TCR alpha chain."""
        score = 0.0
        
        # Length check with stronger discrimination
        alpha_range = self.tcr_patterns['alpha']['length_range']
        optimal_range = self.tcr_patterns['alpha']['optimal_length']
        
        if alpha_range[0] <= length <= alpha_range[1]:
            score += 0.4
            # Bonus for optimal length range
            if optimal_range[0] <= length <= optimal_range[1]:
                score += 0.2
        elif 170 <= length < alpha_range[0]:
            # Slightly short alpha chains are acceptable
            score += 0.3
        elif length > alpha_range[1]:
            # Penalty for being too long (likely beta)
            score -= 0.3
        
        # N-terminal patterns (very distinctive for alpha)
        n_terminal = sequence[:20] if len(sequence) >= 20 else sequence
        n_patterns = self.tcr_patterns['alpha']['n_terminal_patterns']
        for pattern in n_patterns:
            if pattern in n_terminal:
                score += 0.3
                break
        
        # Motif checks
        motifs = self.tcr_patterns['alpha']['conserved_motifs']
        motif_count = sum(1 for motif in motifs if motif in sequence)
        score += motif_count * 0.12
        
        # CDR3 patterns
        cdr3_patterns = self.tcr_patterns['alpha']['cdr3_patterns']
        cdr3_count = sum(1 for pattern in cdr3_patterns if pattern in sequence)
        score += cdr3_count * 0.08
        
        # TCR-specific C-terminal patterns
        c_terminal = sequence[-20:] if len(sequence) >= 20 else sequence
        if any(pattern in c_terminal for pattern in ['FGXGT', 'FGAGT', 'FGQGT']):
            score += 0.2
        
        # Strong preference for shorter chains (alpha characteristic)
        if length < 190:
            score += 0.25
        elif length < 220:
            score += 0.2
        elif length > 250:
            score -= 0.4  # Strong penalty for long chains
        
        return min(max(score, 0.0), 1.0)
    
    def _score_tcr_beta(self, sequence: str, length: int) -> float:
        """Score likelihood of being TCR beta chain."""
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
        
        # CDR3 patterns
        cdr3_patterns = self.tcr_patterns['beta']['cdr3_patterns']
        cdr3_count = sum(1 for pattern in cdr3_patterns if pattern in sequence)
        score += cdr3_count * 0.08
        
        # Strong preference for longer chains (beta characteristic)
        if length > 230 and length < 290:
            score += 0.2
        elif length < 220:
            score -= 0.5  # Strong penalty for short chains
        elif length > 300:
            score -= 0.6  # Very strong penalty for very long chains (likely MHC)
        
        # Penalty for MHC-specific patterns
        mhc_anti_patterns = self.mhc_patterns.get('anti_tcr_patterns', [])
        if any(pattern in sequence for pattern in mhc_anti_patterns):
            score -= 0.8  # Strong penalty for MHC patterns
        
        # Penalty for alpha-specific patterns
        alpha_specific = self.tcr_patterns['alpha']['highly_specific_motifs']
        if any(pattern in sequence for pattern in alpha_specific):
            score -= 0.6  # Penalty for alpha-specific patterns
        
        return min(max(score, 0.0), 1.0)
    
    def _detect_complexes(self, chain_info: Dict, structure) -> List[Dict]:
        """
        Detect multiple pHLA-TCR complexes in the structure.
        
        Args:
            chain_info: Dictionary with chain information
            structure: BioPython structure object
            
        Returns:
            List of dictionaries, each containing chains for one complex
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
        Detect complexes by counting potential chain types.
        """
        # Score all chains for all types
        potential_types = {
            'peptide': [],
            'b2m': [],
            'mhc_heavy': [],
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
                'b2m': self._score_b2m(seq, length, props),
                'mhc_heavy': self._score_mhc(seq, length, props),
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
        
        # Count strong candidates for unique types
        strong_peptides = len([s for _, s in potential_types['peptide'] if s > 0.6])
        strong_b2m = len([s for _, s in potential_types['b2m'] if s > 0.6])
        strong_mhc = len([s for _, s in potential_types['mhc_heavy'] if s > 0.6])
        
        # Estimate number of complexes
        num_complexes = max(strong_peptides, strong_b2m, strong_mhc, 1)
        
        if num_complexes == 1:
            return [chain_info]
        
        # Group chains into complexes
        complexes = []
        assigned_chains = set()
        
        for i in range(num_complexes):
            complex_chains = {}
            
            # Assign one of each unique type per complex, but be more flexible for MHC
            for unique_type in ['peptide', 'b2m', 'mhc_heavy']:
                candidates = potential_types[unique_type]
                score_threshold = 0.3 if unique_type != 'mhc_heavy' else 0.1  # Lower threshold for MHC
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
    
    def _classify_multiple_complexes(self, complexes: List[Dict]) -> Dict[str, str]:
        """
        Classify chains within multiple complexes.
        
        Args:
            complexes: List of dictionaries, each containing chains for one complex
            
        Returns:
            Dictionary mapping all chain IDs to their types
        """
        all_assignments = {}
        
        for i, complex_chains in enumerate(complexes):
            if self.verbose:
                print(f"\nClassifying Complex {i+1} ({len(complex_chains)} chains):")
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
    
    def _validate_assignments(self, assignments: Dict[str, str], 
                            chain_info: Dict) -> Dict[str, str]:
        """
        Validate and adjust chain assignments to ensure consistency.
        
        Args:
            assignments: Initial assignments
            chain_info: Chain information
            
        Returns:
            Validated assignments
        """
        # Count each type
        type_counts = defaultdict(int)
        for chain_type in assignments.values():
            type_counts[chain_type] += 1
        
        # Ensure we don't have duplicates of unique chains
        unique_types = ['mhc_heavy', 'b2m', 'peptide']
        for unique_type in unique_types:
            if type_counts[unique_type] > 1:
                # Keep the best candidate, reassign others
                candidates = [(cid, info) for cid, info in chain_info.items() 
                             if assignments[cid] == unique_type]
                
                # Keep the best candidate based on scoring
                best_candidate = None
                best_score = -1
                
                for chain_id, info in candidates:
                    if unique_type == 'peptide':
                        score = self._score_peptide(info['sequence'], info['length'], info['properties'])
                    elif unique_type == 'b2m':
                        score = self._score_b2m(info['sequence'], info['length'], info['properties'])
                    elif unique_type == 'mhc_heavy':
                        score = self._score_mhc(info['sequence'], info['length'], info['properties'])
                    else:
                        score = 0
                    
                    if score > best_score:
                        best_score = score
                        best_candidate = chain_id
                
                # Reassign all but the best candidate
                for chain_id, _ in candidates:
                    if chain_id != best_candidate:
                        assignments[chain_id] = 'unknown'
        
        return assignments
    
    def save_analysis_report(self, assignments: Dict[str, str], 
                           output_file: str, chain_info: Dict = None):
        """
        Save analysis report to file.
        
        Args:
            assignments: Chain assignments
            output_file: Output file path
            chain_info: Additional chain information
        """
        with open(output_file, 'w') as f:
            f.write("pHLA-TCR Complex Analysis Report\n")
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


def main():
    """Command-line interface for the analyzer."""
    import argparse
    
    parser = argparse.ArgumentParser(
        description="Analyze pHLA-TCR complex PDB structures"
    )
    parser.add_argument("pdb_file", help="Input PDB file")
    parser.add_argument("-o", "--output", help="Output report file")
    parser.add_argument("-v", "--verbose", action="store_true", 
                       help="Verbose output")
    
    args = parser.parse_args()
    
    # Initialize analyzer
    analyzer = pMHCITCRAnalyzer(verbose=args.verbose)
    
    try:
        # Analyze the PDB file
        assignments = analyzer.analyze_pdb(args.pdb_file)
        
        # Print results in alphabetical order
        print("\nFinal Chain Assignments:")
        print("-" * 30)
        for chain_id in sorted(assignments.keys()):
            print(f"Chain {chain_id}: {assignments[chain_id]}")
        
        # Save report if requested
        if args.output:
            analyzer.save_analysis_report(assignments, args.output)
            print(f"\nReport saved to: {args.output}")
            
    except Exception as e:
        print(f"Error: {e}")
        sys.exit(1)


if __name__ == "__main__":
    main()