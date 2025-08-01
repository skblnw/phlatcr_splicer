#!/usr/bin/env python3
"""
Test script for pHLA-TCR Complex Structure Analyzer

This script provides test cases and validation for the analyzer.
"""

import unittest
import tempfile
import os
from phlatcr_analyzer import pHLATCRAnalyzer


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
        # Mock sequence for TCR alpha
        alpha_seq = "EVQLVESGGGLVQPGGSLRLSCAASGFTFSSYGMHWVRQAPGKGLEWVAVISYDGSNKYYADSVKGRFTISRDNSKNTVYLQMNSLRAEDTAVYYCAKFGXGTQLTVSS"
        alpha_score = self.analyzer._score_tcr_alpha(alpha_seq, len(alpha_seq))
        self.assertGreater(alpha_score, 0)
        
        # Mock sequence for TCR beta
        beta_seq = "MGIGVTQTPKFQVLKTGQSMTLQCAQDMNHEYMSWYRQDPGMGLRLIHYSVGAGTDKGEVPNPFDSAQGVKLISEDGYDYICHFFGTGTRLTVLEDLKNV"
        beta_score = self.analyzer._score_tcr_beta(beta_seq, len(beta_seq))
        self.assertGreater(beta_score, 0)
    
    def test_chain_patterns(self):
        """Test pattern loading."""
        self.assertIn('length_range', self.analyzer.mhc_patterns)
        self.assertIn('alpha', self.analyzer.tcr_patterns)
        self.assertIn('beta', self.analyzer.tcr_patterns)
        self.assertIn('conserved_motifs', self.analyzer.b2m_patterns)
    
    def create_mock_pdb(self) -> str:
        """Create a minimal mock PDB file for testing."""
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
        # Create temporary file
        temp_file = tempfile.NamedTemporaryFile(mode='w', suffix='.pdb', delete=False)
        temp_file.write(pdb_content)
        temp_file.close()
        return temp_file.name
    
    def test_pdb_parsing(self):
        """Test basic PDB parsing functionality."""
        mock_pdb = self.create_mock_pdb()
        
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
    
    def test_file_not_found(self):
        """Test handling of non-existent files."""
        with self.assertRaises(FileNotFoundError):
            self.analyzer.analyze_pdb("nonexistent_file.pdb")


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