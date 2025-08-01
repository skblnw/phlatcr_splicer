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
    from Bio.PDB import PDBParser, Structure, Model, Chain, Residue
    from Bio.PDB.PDBIO import PDBIO
    from Bio.SeqUtils import ProtParam
    from Bio.Seq import Seq
except ImportError:
    print("Error: BioPython is required. Install with: pip install biopython")
    sys.exit(1)


class pHLATCRAnalyzer:
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
            'length_range': (270, 380),  # Typical length range
            'conserved_motifs': [
                'GSHSMRY',  # Alpha1 domain motif
                'WSDRVI',   # Alpha2 domain motif
                'FKAFLKQ',  # Alpha3 domain motif
            ],
            'domain_boundaries': {
                'alpha1': (1, 90),
                'alpha2': (91, 180),
                'alpha3': (181, 276),
            }
        }
    
    def _load_tcr_patterns(self) -> Dict:
        """Load TCR alpha and beta chain patterns."""
        return {
            'alpha': {
                'length_range': (200, 250),
                'conserved_motifs': [
                    'FGXGT',    # CDR3 C-terminal
                    'WYQQKP',   # Framework regions
                ],
                'v_gene_patterns': ['TRAV', 'TCRAV'],
                'j_gene_patterns': ['TRAJ', 'TCRAJ'],
            },
            'beta': {
                'length_range': (240, 290),
                'conserved_motifs': [
                    'FGGGT',    # CDR3 C-terminal
                    'WYQQKP',   # Framework regions
                ],
                'v_gene_patterns': ['TRBV', 'TCRBV'],
                'j_gene_patterns': ['TRBJ', 'TCRBJ'],
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
        
        # Classify chains
        chain_assignments = self._classify_chains(chain_info)
        
        if self.verbose:
            print("\nChain assignments:")
            for chain_id, chain_type in chain_assignments.items():
                print(f"  Chain {chain_id}: {chain_type}")
        
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
        
        for chain_id, info in sorted_chains:
            sequence = info['sequence']
            length = info['length']
            properties = info['properties']
            
            # Classification logic
            chain_type = self._determine_chain_type(
                sequence, length, properties, used_types
            )
            
            assignments[chain_id] = chain_type
            used_types.add(chain_type)
        
        # Post-process to ensure we have reasonable assignments
        assignments = self._validate_assignments(assignments, chain_info)
        
        return assignments
    
    def _determine_chain_type(self, sequence: str, length: int, 
                            properties: Dict, used_types: Set) -> str:
        """
        Determine the type of a single chain.
        
        Args:
            sequence: Amino acid sequence
            length: Sequence length
            properties: Chain properties
            used_types: Already assigned chain types
            
        Returns:
            Predicted chain type
        """
        # β2-microglobulin is usually shortest and has distinctive pattern
        if (self.b2m_patterns['length_range'][0] <= length <= 
            self.b2m_patterns['length_range'][1] and 
            'b2m' not in used_types):
            
            # Check for β2m conserved motifs
            b2m_score = sum(1 for motif in self.b2m_patterns['conserved_motifs'] 
                           if motif in sequence)
            if b2m_score >= 1:  # At least one conserved motif
                return 'b2m'
        
        # Peptide antigen (usually shortest after β2m)
        if length < 20 and 'peptide' not in used_types:
            return 'peptide'
        
        # MHC heavy chain (typically longest)
        if (self.mhc_patterns['length_range'][0] <= length <= 
            self.mhc_patterns['length_range'][1] and 
            'mhc_heavy' not in used_types):
            
            # Check for MHC conserved motifs
            mhc_score = sum(1 for motif in self.mhc_patterns['conserved_motifs'] 
                           if motif in sequence)
            if mhc_score >= 1:
                return 'mhc_heavy'
        
        # TCR chains
        if 'tcr_alpha' not in used_types:
            alpha_score = self._score_tcr_alpha(sequence, length)
            if alpha_score > 0.5:
                return 'tcr_alpha'
        
        if 'tcr_beta' not in used_types:
            beta_score = self._score_tcr_beta(sequence, length)
            if beta_score > 0.5:
                return 'tcr_beta'
        
        # Default classification based on length and what's missing
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
        
        # Length check
        alpha_range = self.tcr_patterns['alpha']['length_range']
        if alpha_range[0] <= length <= alpha_range[1]:
            score += 0.3
        
        # Motif checks
        motifs = self.tcr_patterns['alpha']['conserved_motifs']
        for motif in motifs:
            if motif in sequence:
                score += 0.2
        
        # Additional TCR-specific patterns
        if 'FGXGT' in sequence or 'FGGGT' in sequence:
            score += 0.3
        
        return min(score, 1.0)
    
    def _score_tcr_beta(self, sequence: str, length: int) -> float:
        """Score likelihood of being TCR beta chain."""
        score = 0.0
        
        # Length check
        beta_range = self.tcr_patterns['beta']['length_range']
        if beta_range[0] <= length <= beta_range[1]:
            score += 0.3
        
        # Motif checks
        motifs = self.tcr_patterns['beta']['conserved_motifs']
        for motif in motifs:
            if motif in sequence:
                score += 0.2
        
        # TCR beta tends to be longer than alpha
        if length > 250:
            score += 0.2
        
        return min(score, 1.0)
    
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
        unique_types = ['mhc_heavy', 'b2m']
        for unique_type in unique_types:
            if type_counts[unique_type] > 1:
                # Keep the best candidate, reassign others
                candidates = [(cid, info) for cid, info in chain_info.items() 
                             if assignments[cid] == unique_type]
                
                # Keep the best candidate (could implement better scoring here)
                best_candidate = max(candidates, key=lambda x: x[1]['length'])
                
                for chain_id, _ in candidates:
                    if chain_id != best_candidate[0]:
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
            
            for chain_id, chain_type in assignments.items():
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
    analyzer = pHLATCRAnalyzer(verbose=args.verbose)
    
    try:
        # Analyze the PDB file
        assignments = analyzer.analyze_pdb(args.pdb_file)
        
        # Print results
        print("\nFinal Chain Assignments:")
        print("-" * 30)
        for chain_id, chain_type in assignments.items():
            print(f"Chain {chain_id}: {chain_type}")
        
        # Save report if requested
        if args.output:
            analyzer.save_analysis_report(assignments, args.output)
            print(f"\nReport saved to: {args.output}")
            
    except Exception as e:
        print(f"Error: {e}")
        sys.exit(1)


if __name__ == "__main__":
    main()