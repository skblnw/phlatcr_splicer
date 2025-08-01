#!/usr/bin/env python3
"""
Test script for pHLA-TCR Complex Structure Analyzer

This script provides test cases and validation for the analyzer.
"""

import unittest
import tempfile
import os
import sys
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..'))
from phlatcr_splicer import pHLATCRAnalyzer


class TestpHLATCRAnalyzer(unittest.TestCase):
    """Test cases for the pHLA-TCR analyzer."""
    
    def setUp(self):
        """Set up test fixtures."""
        self.analyzer = pHLATCRAnalyzer(verbose=False)
    
    def test_three_to_one_conversion(self):
        """Test amino acid 3-letter to 1-letter conversion."""
        self.assertEqual(self.analyzer._three_to_one('ALA'), 'A')
        self.assertEqual(self.analyzer._three_to_one('GLY'), 'G')
        self.assertEqual(self.analyzer._three_to_one('TRP'), 'W')
        self.assertEqual(self.analyzer._three_to_one('XXX'), '')
    
    def test_tcr_scoring(self):
        """Test TCR scoring functions."""
        # Mock sequence for TCR alpha (based on real patterns)
        alpha_seq = "QLLEQSPQFLSIQEGENLTVYCNSSSVFSSLQWYRQEPGEGPVLLVTVVT" + "A" * 150  # 200 residues
        alpha_score = self.analyzer._score_tcr_alpha(alpha_seq, len(alpha_seq))
        self.assertGreater(alpha_score, 0.5)
        
        # Mock sequence for TCR beta (based on real patterns)  
        beta_seq = "GITQSPKYLFRKEGQNVTLSCEQNLNHDAMYWYRQDPGQGLRLIYYSQIV" + "B" * 190  # 240 residues
        beta_score = self.analyzer._score_tcr_beta(beta_seq, len(beta_seq))
        self.assertGreater(beta_score, 0.5)
    
    def test_chain_patterns(self):
        """Test pattern loading."""
        self.assertIn('length_range', self.analyzer.mhc_patterns)
        self.assertIn('alpha', self.analyzer.tcr_patterns)
        self.assertIn('beta', self.analyzer.tcr_patterns)
        self.assertIn('conserved_motifs', self.analyzer.b2m_patterns)
    
    def create_mock_pdb(self, complex_type="simple") -> str:
        """Create mock PDB files for testing different scenarios."""
        
        if complex_type == "simple":
            # Simple 2-chain complex
            pdb_content = """HEADER    IMMUNE SYSTEM/IMMUNE SYSTEM             01-JAN-20   TEST    
ATOM      1  N   MET A   1      20.154  16.000  12.000  1.00 10.00           N  
ATOM      2  CA  MET A   1      19.000  16.500  12.500  1.00 10.00           C  
ATOM      3  C   MET A   1      18.000  15.500  13.000  1.00 10.00           C  
ATOM      4  O   MET A   1      17.500  15.000  12.000  1.00 10.00           O  
ATOM      5  N   GLY B   1      25.000  20.000  15.000  1.00 10.00           N  
ATOM      6  CA  GLY B   1      24.500  20.500  15.500  1.00 10.00           C  
ATOM      7  C   GLY B   1      23.500  19.500  16.000  1.00 10.00           C  
ATOM      8  O   GLY B   1      23.000  19.000  15.000  1.00 10.00           O  
END
"""
        elif complex_type == "phla_tcr":
            # Mock pHLA-TCR complex with realistic chain lengths
            pdb_content = """HEADER    IMMUNE SYSTEM/IMMUNE SYSTEM             01-JAN-20   MOCK    
TITLE     MOCK PHLA-TCR COMPLEX
"""
            # Add peptide chain (chain P) - 9 residues
            residue_num = 1
            atom_num = 1
            for i, aa in enumerate(['TYR', 'LEU', 'GLN', 'PRO', 'ARG', 'TRP', 'TYR', 'PHE', 'VAL']):
                x, y, z = 20.0 + i, 15.0, 10.0
                pdb_content += f"ATOM  {atom_num:5d}  N   {aa} P{i+1:4d}    {x:8.3f}{y:8.3f}{z:8.3f}  1.00 20.00           N  \n"
                atom_num += 1
                pdb_content += f"ATOM  {atom_num:5d}  CA  {aa} P{i+1:4d}    {x+0.5:8.3f}{y+0.5:8.3f}{z+0.5:8.3f}  1.00 20.00           C  \n"
                atom_num += 1
            
            # Add β2m chain (chain B) - 99 residues (simplified)
            for i in range(99):
                aa = 'ALA' if i % 2 == 0 else 'GLY'
                x, y, z = 30.0 + i*0.1, 20.0, 15.0
                pdb_content += f"ATOM  {atom_num:5d}  N   {aa} B{i+1:4d}    {x:8.3f}{y:8.3f}{z:8.3f}  1.00 20.00           N  \n"
                atom_num += 1
                pdb_content += f"ATOM  {atom_num:5d}  CA  {aa} B{i+1:4d}    {x+0.5:8.3f}{y+0.5:8.3f}{z+0.5:8.3f}  1.00 20.00           C  \n"
                atom_num += 1
            
            # Add MHC heavy chain (chain H) - 276 residues (simplified)
            for i in range(276):
                aa = ['ALA', 'VAL', 'LEU', 'ILE'][i % 4]
                x, y, z = 40.0 + i*0.1, 25.0, 20.0
                pdb_content += f"ATOM  {atom_num:5d}  N   {aa} H{i+1:4d}    {x:8.3f}{y:8.3f}{z:8.3f}  1.00 20.00           N  \n"
                atom_num += 1
                pdb_content += f"ATOM  {atom_num:5d}  CA  {aa} H{i+1:4d}    {x+0.5:8.3f}{y+0.5:8.3f}{z+0.5:8.3f}  1.00 20.00           C  \n"
                atom_num += 1
            
            # Add TCR alpha chain (chain A) - 220 residues
            for i in range(220):
                aa = ['SER', 'THR', 'ASN', 'GLN'][i % 4]
                x, y, z = 50.0 + i*0.1, 30.0, 25.0
                pdb_content += f"ATOM  {atom_num:5d}  N   {aa} A{i+1:4d}    {x:8.3f}{y:8.3f}{z:8.3f}  1.00 20.00           N  \n"
                atom_num += 1
                pdb_content += f"ATOM  {atom_num:5d}  CA  {aa} A{i+1:4d}    {x+0.5:8.3f}{y+0.5:8.3f}{z+0.5:8.3f}  1.00 20.00           C  \n"
                atom_num += 1
            
            # Add TCR beta chain (chain T) - 245 residues
            for i in range(245):
                aa = ['GLU', 'ASP', 'LYS', 'ARG'][i % 4]
                x, y, z = 60.0 + i*0.1, 35.0, 30.0
                pdb_content += f"ATOM  {atom_num:5d}  N   {aa} T{i+1:4d}    {x:8.3f}{y:8.3f}{z:8.3f}  1.00 20.00           N  \n"
                atom_num += 1
                pdb_content += f"ATOM  {atom_num:5d}  CA  {aa} T{i+1:4d}    {x+0.5:8.3f}{y+0.5:8.3f}{z+0.5:8.3f}  1.00 20.00           C  \n"
                atom_num += 1
            
            pdb_content += "END\n"
        
        # Create temporary file
        temp_file = tempfile.NamedTemporaryFile(mode='w', suffix='.pdb', delete=False)
        temp_file.write(pdb_content)
        temp_file.close()
        return temp_file.name
    
    def test_pdb_parsing(self):
        """Test basic PDB parsing functionality."""
        mock_pdb = self.create_mock_pdb("simple")
        
        try:
            # Test that the analyzer can parse the file without errors
            result = self.analyzer.analyze_pdb(mock_pdb)
            self.assertIsInstance(result, dict)
            
            # Should have two chains (A and B)
            self.assertEqual(len(result), 2)
            self.assertIn('A', result)
            self.assertIn('B', result)
            
        finally:
            # Clean up
            os.unlink(mock_pdb)
    
    def test_phla_tcr_complex(self):
        """Test analysis of a mock pHLA-TCR complex."""
        mock_pdb = self.create_mock_pdb("phla_tcr")
        
        try:
            result = self.analyzer.analyze_pdb(mock_pdb)
            self.assertIsInstance(result, dict)
            
            # Should have 5 chains
            self.assertEqual(len(result), 5)
            
            # Check that we have all expected chain types
            chain_types = set(result.values())
            expected_types = {'peptide', 'b2m', 'mhc_heavy', 'tcr_alpha', 'tcr_beta'}
            
            # At least 3 of the 5 expected types should be identified
            overlap = len(chain_types.intersection(expected_types))
            self.assertGreaterEqual(overlap, 3)
            
            if self.analyzer.verbose:
                print(f"Mock pHLA-TCR analysis result: {result}")
            
        finally:
            os.unlink(mock_pdb)
    
    def test_file_not_found(self):
        """Test handling of non-existent files."""
        with self.assertRaises(FileNotFoundError):
            self.analyzer.analyze_pdb("nonexistent_file.pdb")
    
    def test_scoring_functions(self):
        """Test individual scoring functions."""
        # Test peptide scoring
        peptide_seq = "YLQPRTWY"  # 8-mer peptide
        peptide_score = self.analyzer._score_peptide(peptide_seq, len(peptide_seq), {})
        self.assertGreater(peptide_score, 0.5)
        
        # Test long sequence (shouldn't be peptide)
        long_seq = "A" * 100
        long_score = self.analyzer._score_peptide(long_seq, len(long_seq), {})
        self.assertLess(long_score, 0.5)
        
        # Test β2m scoring
        b2m_seq = "A" * 99  # Right length
        b2m_score = self.analyzer._score_b2m(b2m_seq, len(b2m_seq), {'molecular_weight': 12000})
        self.assertGreater(b2m_score, 0.5)
        
        # Test MHC scoring
        mhc_seq = "A" * 320  # Right length for MHC
        mhc_score = self.analyzer._score_mhc(mhc_seq, len(mhc_seq), {})
        self.assertGreater(mhc_score, 0.3)
    
    def test_chain_validation(self):
        """Test chain assignment validation."""
        # Create mock chain info
        chain_info = {
            'A': {'length': 10, 'sequence': 'YLQPRTWYDF', 'properties': {}},
            'B': {'length': 99, 'sequence': 'A' * 99, 'properties': {'molecular_weight': 12000}},
            'C': {'length': 320, 'sequence': 'G' * 320, 'properties': {}},
        }
        
        # Mock initial assignments (with duplicates)
        initial_assignments = {'A': 'peptide', 'B': 'peptide', 'C': 'mhc_heavy'}
        
        # Validate
        validated = self.analyzer._validate_assignments(initial_assignments, chain_info)
        
        # Should resolve duplicates
        peptide_count = sum(1 for v in validated.values() if v == 'peptide')
        self.assertLessEqual(peptide_count, 1)


def run_example_analysis():
    """Run an example analysis with mock data."""
    print("Running example analysis...")
    
    analyzer = pHLATCRAnalyzer(verbose=True)
    
    # Create test data
    test = TestpHLATCRAnalyzer()
    test.setUp()
    mock_pdb = test.create_mock_pdb()
    
    try:
        print(f"\nAnalyzing mock PDB file: {mock_pdb}")
        result = analyzer.analyze_pdb(mock_pdb)
        
        print("\nResults:")
        for chain_id, chain_type in result.items():
            print(f"  Chain {chain_id}: {chain_type}")
        
        # Save report
        report_file = "example_analysis_report.txt"
        analyzer.save_analysis_report(result, report_file)
        print(f"\nReport saved to: {report_file}")
        
    finally:
        os.unlink(mock_pdb)


if __name__ == "__main__":
    print("pHLA-TCR Analyzer Test Suite")
    print("=" * 40)
    
    # Run unit tests
    print("\n1. Running unit tests...")
    unittest.main(argv=[''], exit=False, verbosity=2)
    
    # Run example analysis
    print("\n2. Running example analysis...")
    run_example_analysis()
    
    print("\nTest suite completed!")