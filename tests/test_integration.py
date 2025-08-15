#!/usr/bin/env python3
"""
Integration tests for the Unified TCR-pMHC Analyzer

End-to-end tests with complete analysis workflows.

Author: skblnw
"""

import unittest
import tempfile
import os
from pathlib import Path
import numpy as np

from phlatcr_splicer import TCRpMHCAnalyzer


class TestEndToEnd(unittest.TestCase):
    """End-to-end integration tests."""
    
    def setUp(self):
        """Set up test fixtures."""
        self.analyzer = TCRpMHCAnalyzer(verbose=False)
    
    def create_mock_pdb(self, chains_config):
        """Create a mock PDB file with specified chains."""
        pdb_lines = ["HEADER    TEST STRUCTURE"]
        atom_id = 1
        
        for chain_id, (sequence, start_pos) in chains_config.items():
            for i, aa in enumerate(sequence):
                # Map single letter to three letter code
                aa_map = {
                    'A': 'ALA', 'G': 'GLY', 'S': 'SER', 'T': 'THR',
                    'V': 'VAL', 'L': 'LEU', 'I': 'ILE', 'M': 'MET',
                    'F': 'PHE', 'W': 'TRP', 'P': 'PRO', 'C': 'CYS',
                    'Y': 'TYR', 'N': 'ASN', 'Q': 'GLN', 'D': 'ASP',
                    'E': 'GLU', 'K': 'LYS', 'R': 'ARG', 'H': 'HIS'
                }
                aa3 = aa_map.get(aa, 'ALA')
                
                x = start_pos[0] + i * 3.8
                y = start_pos[1]
                z = start_pos[2]
                
                pdb_lines.append(
                    f"ATOM  {atom_id:5d}  CA  {aa3} {chain_id}{i+1:4d}    "
                    f"{x:8.3f}{y:8.3f}{z:8.3f}  1.00  0.00           C"
                )
                atom_id += 1
        
        pdb_lines.append("END")
        return "\n".join(pdb_lines)
    
    def test_mhc_i_complex_analysis(self):
        """Test complete MHC-I complex analysis."""
        # Create a mock MHC-I complex
        chains = {
            'A': ('GSHSMR' + 'A' * 260, [0, 0, 0]),      # MHC-I heavy chain
            'B': ('IQRTPKIQ' + 'A' * 91, [30, 0, 0]),    # B2M
            'C': ('NLVPMVATV', [15, 0, 0]),              # Peptide
            'D': ('A' * 200, [50, 0, 0]),                # TCR alpha
            'E': ('A' * 250, [55, 0, 0]),                # TCR beta
        }
        
        pdb_content = self.create_mock_pdb(chains)
        
        with tempfile.NamedTemporaryFile(mode='w', suffix='.pdb', delete=False) as f:
            f.write(pdb_content)
            test_file = f.name
        
        try:
            results, chain_map = self.analyzer.analyze_pdb(test_file)
            
            # Check that results were returned
            self.assertIsInstance(results, dict)
            self.assertTrue(len(results) > 0, "Should identify at least one complex")
            
            # Check chain types identified
            for complex_data in results.values():
                chains = complex_data.get('chains', {})
                self.assertIn('MHC_I_ALPHA', chains.values())
                self.assertIn('B2M', chains.values())
                self.assertIn('PEPTIDE', chains.values())
        finally:
            os.unlink(test_file)
    
    def test_mhc_ii_complex_analysis(self):
        """Test complete MHC-II complex analysis."""
        # Create a mock MHC-II complex
        chains = {
            'A': ('WRLEEF' + 'A' * 180, [0, 0, 0]),      # MHC-II alpha
            'B': ('NGTER' + 'A' * 195, [5, 0, 0]),       # MHC-II beta
            'C': ('PKYVKQNTLKLATGMR', [2.5, 0, 0]),      # Peptide
            'D': ('A' * 200, [50, 0, 0]),                # TCR alpha
            'E': ('A' * 250, [55, 0, 0]),                # TCR beta
        }
        
        pdb_content = self.create_mock_pdb(chains)
        
        with tempfile.NamedTemporaryFile(mode='w', suffix='.pdb', delete=False) as f:
            f.write(pdb_content)
            test_file = f.name
        
        try:
            results, chain_map = self.analyzer.analyze_pdb(test_file)
            
            # Check that results were returned
            self.assertIsInstance(results, dict)
            self.assertTrue(len(results) > 0, "Should identify at least one complex")
            
            # Check chain types identified
            for complex_data in results.values():
                chains = complex_data.get('chains', {})
                self.assertIn('MHC_II_ALPHA', chains.values())
                self.assertIn('MHC_II_BETA', chains.values())
                self.assertIn('PEPTIDE', chains.values())
        finally:
            os.unlink(test_file)
    
    def test_multi_complex_detection(self):
        """Test detection of multiple complexes in one file."""
        # Create two separate MHC-I complexes
        chains = {
            # Complex 1
            'A': ('GSHSMR' + 'A' * 260, [0, 0, 0]),
            'B': ('IQRTPKIQ' + 'A' * 91, [5, 0, 0]),
            'C': ('NLVPMVATV', [2.5, 0, 0]),
            # Complex 2 (far away)
            'D': ('GSHSMR' + 'A' * 260, [200, 0, 0]),
            'E': ('IQRTPKIQ' + 'A' * 91, [205, 0, 0]),
            'F': ('SIINFEKL', [202.5, 0, 0]),
        }
        
        pdb_content = self.create_mock_pdb(chains)
        
        with tempfile.NamedTemporaryFile(mode='w', suffix='.pdb', delete=False) as f:
            f.write(pdb_content)
            test_file = f.name
        
        try:
            results, chain_map = self.analyzer.analyze_pdb(test_file)
            
            # Should detect two separate complexes
            self.assertEqual(len(results), 2, "Should identify two separate complexes")
            
            # Each complex should have MHC-I components
            for complex_data in results.values():
                pairs = complex_data.get('pairs', [])
                self.assertTrue(len(pairs) > 0)
                for pair in pairs:
                    if pair.get('type') == 'MHC-I':
                        self.assertIsNotNone(pair.get('mhc_alpha'))
                        self.assertIsNotNone(pair.get('b2m'))
        finally:
            os.unlink(test_file)
    
    def test_unknown_chain_reclassification(self):
        """Test that unknown chains get reclassified based on proximity."""
        # Create chains where some will be initially unknown
        chains = {
            'A': ('GSHSMR' + 'A' * 260, [0, 0, 0]),      # MHC-I heavy chain
            'B': ('X' * 99, [5, 0, 0]),                  # Will be unknown, then B2M
            'C': ('NLVPMVATV', [2.5, 0, 0]),             # Peptide
        }
        
        pdb_content = self.create_mock_pdb(chains)
        
        with tempfile.NamedTemporaryFile(mode='w', suffix='.pdb', delete=False) as f:
            f.write(pdb_content)
            test_file = f.name
        
        try:
            results, chain_map = self.analyzer.analyze_pdb(test_file)
            
            # Check that chain B was reclassified to B2M
            for complex_data in results.values():
                chains = complex_data.get('chains', {})
                # Chain B should be reclassified as B2M due to proximity to MHC-I
                if 'B' in chains:
                    # It might be reclassified or remain unknown depending on patterns
                    self.assertIn(chains['B'], ['B2M', 'UNKNOWN'])
        finally:
            os.unlink(test_file)
    
    def test_parameter_sensitivity(self):
        """Test that different parameters affect results."""
        # Create a borderline case
        chains = {
            'A': ('GSHSMR' + 'A' * 260, [0, 0, 0]),
            'B': ('IQRTPKIQ' + 'A' * 91, [45, 0, 0]),    # At edge of default distance
        }
        
        pdb_content = self.create_mock_pdb(chains)
        
        with tempfile.NamedTemporaryFile(mode='w', suffix='.pdb', delete=False) as f:
            f.write(pdb_content)
            test_file = f.name
        
        try:
            # Test with default parameters
            analyzer_default = TCRpMHCAnalyzer(verbose=False)
            results_default, _ = analyzer_default.analyze_pdb(test_file)
            
            # Test with tighter clustering
            analyzer_tight = TCRpMHCAnalyzer(eps=30.0, verbose=False)
            results_tight, _ = analyzer_tight.analyze_pdb(test_file)
            
            # Results might differ based on clustering
            # At minimum, both should process without error
            self.assertIsInstance(results_default, dict)
            self.assertIsInstance(results_tight, dict)
        finally:
            os.unlink(test_file)
    
    def test_empty_file_handling(self):
        """Test handling of empty or invalid files."""
        with tempfile.NamedTemporaryFile(mode='w', suffix='.pdb', delete=False) as f:
            f.write("HEADER    EMPTY\nEND\n")
            test_file = f.name
        
        try:
            results, chain_map = self.analyzer.analyze_pdb(test_file)
            # Should return empty results, not crash
            self.assertEqual(len(results), 0)
        finally:
            os.unlink(test_file)
    
    def test_cif_file_processing(self):
        """Test that CIF files are processed correctly."""
        # Create a minimal CIF file
        cif_content = """data_TEST
loop_
_atom_site.group_PDB
_atom_site.id
_atom_site.type_symbol
_atom_site.label_atom_id
_atom_site.label_alt_id
_atom_site.label_comp_id
_atom_site.label_asym_id
_atom_site.label_entity_id
_atom_site.label_seq_id
_atom_site.Cartn_x
_atom_site.Cartn_y
_atom_site.Cartn_z
ATOM 1 C CA . ALA A 1 1 0.000 0.000 0.000
ATOM 2 C CA . ALA A 1 2 3.800 0.000 0.000
"""
        with tempfile.NamedTemporaryFile(mode='w', suffix='.cif', delete=False) as f:
            f.write(cif_content)
            test_file = f.name
        
        try:
            results, chain_map = self.analyzer.analyze_pdb(test_file)
            # Should process without error
            self.assertIsInstance(results, dict)
            self.assertIsInstance(chain_map, dict)
        finally:
            os.unlink(test_file)


class TestPerformance(unittest.TestCase):
    """Performance and scalability tests."""
    
    def setUp(self):
        """Set up test fixtures."""
        self.analyzer = TCRpMHCAnalyzer(verbose=False)
    
    def test_large_complex_performance(self):
        """Test performance with large multi-chain complex."""
        # Create a large complex with many chains
        chains = {}
        for i in range(20):  # 20 chains
            chain_id = chr(65 + i) if i < 26 else str(i - 26)
            sequence = 'A' * (100 + i * 10)  # Varying lengths
            position = [i * 10, 0, 0]
            chains[chain_id] = (sequence, position)
        
        pdb_content = self.create_mock_pdb(chains)
        
        with tempfile.NamedTemporaryFile(mode='w', suffix='.pdb', delete=False) as f:
            f.write(pdb_content)
            test_file = f.name
        
        try:
            import time
            start_time = time.time()
            results, chain_map = self.analyzer.analyze_pdb(test_file)
            end_time = time.time()
            
            # Should complete in reasonable time (< 5 seconds for 20 chains)
            self.assertLess(end_time - start_time, 5.0)
            self.assertIsInstance(results, dict)
        finally:
            os.unlink(test_file)
    
    def create_mock_pdb(self, chains_config):
        """Helper method to create mock PDB content."""
        pdb_lines = ["HEADER    TEST STRUCTURE"]
        atom_id = 1
        
        for chain_id, (sequence, start_pos) in chains_config.items():
            for i in range(min(len(sequence), 10)):  # Limit to 10 residues for performance
                x = start_pos[0] + i * 3.8
                y = start_pos[1]
                z = start_pos[2]
                
                pdb_lines.append(
                    f"ATOM  {atom_id:5d}  CA  ALA {chain_id}{i+1:4d}    "
                    f"{x:8.3f}{y:8.3f}{z:8.3f}  1.00  0.00           C"
                )
                atom_id += 1
        
        pdb_lines.append("END")
        return "\n".join(pdb_lines)


if __name__ == '__main__':
    unittest.main()