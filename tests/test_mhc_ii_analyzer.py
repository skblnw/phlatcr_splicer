"""
Unit tests for MHC-II Complex Analyzer

Tests the pMHCIITCRAnalyzer class functionality including:
- MHC-II alpha and beta chain identification
- Longer peptide recognition (12-25 residues)
- TCR alpha/beta chain classification
- Multi-complex detection
- Pattern recognition and scoring

Author: AI Assistant
"""

import unittest
import tempfile
import os
from phlatcr_splicer import pMHCIITCRAnalyzer


class TestpMHCIITCRAnalyzer(unittest.TestCase):
    """Test cases for MHC-II analyzer."""
    
    def setUp(self):
        """Set up test fixtures."""
        self.analyzer = pMHCIITCRAnalyzer(verbose=False)
    
    def test_analyzer_initialization(self):
        """Test analyzer initialization."""
        self.assertIsInstance(self.analyzer, pMHCIITCRAnalyzer)
        self.assertIn('alpha', self.analyzer.mhc_ii_patterns)
        self.assertIn('beta', self.analyzer.mhc_ii_patterns)
        self.assertIn('alpha', self.analyzer.tcr_patterns)
        self.assertIn('beta', self.analyzer.tcr_patterns)
    
    def test_three_to_one_conversion(self):
        """Test amino acid 3-letter to 1-letter conversion."""
        self.assertEqual(self.analyzer._three_to_one('ALA'), 'A')
        self.assertEqual(self.analyzer._three_to_one('VAL'), 'V')
        self.assertEqual(self.analyzer._three_to_one('XXX'), 'X')  # Unknown
    
    def test_peptide_scoring(self):
        """Test MHC-II peptide scoring (longer peptides)."""
        # Short MHC-II peptide (15 residues)
        short_peptide = "PKYVKQNTLKLAT"
        short_score = self.analyzer._score_peptide(short_peptide, len(short_peptide), 
                                                  {'molecular_weight': 1500})
        self.assertGreater(short_score, 0.7)
        
        # Long MHC-II peptide (20 residues)
        long_peptide = "PKYVKQNTLKLATGMRNVPEK"
        long_score = self.analyzer._score_peptide(long_peptide, len(long_peptide), 
                                                 {'molecular_weight': 2200})
        self.assertGreater(long_score, 0.7)
        
        # Too long for peptide
        too_long = "A" * 50
        too_long_score = self.analyzer._score_peptide(too_long, len(too_long), 
                                                     {'molecular_weight': 5000})
        self.assertLess(too_long_score, 0.3)
    
    def test_mhc_ii_alpha_scoring(self):
        """Test MHC-II alpha chain scoring."""
        # Mock sequence with MHC-II alpha patterns
        alpha_seq = "MGSGWVPPLVVSTQLLILHLHPEVQGEDRTQDVTMENEVKNVNLLTLRNGYYNQSAGSHTLQWMHGCELFPFLGPQGEGKASVLVNQNKEMQLFPGTLSGYPGAIVHIVQRLQNTLQKFKDGSLHSLTETLGDFAAHVLRFDSDK" + "A" * 50
        alpha_score = self.analyzer._score_mhc_ii_alpha(alpha_seq, len(alpha_seq))
        self.assertGreater(alpha_score, 0.5)
    
    def test_mhc_ii_beta_scoring(self):
        """Test MHC-II beta chain scoring."""
        # Mock sequence with MHC-II beta patterns  
        beta_seq = "MSWKKALRIPRLQWLLLLTTLLLTTESTQGSQKPTEAGEPEVVTVPVTYQRPVSDLPVQCDIRQTHPGVLRVPHYNQREETEQGYLKEFLNGTERVRFLNGQKRQPKADAIQKLIPKVKRGCHSLHNTLGDFAEKGYVVVTKR" + "A" * 50
        beta_score = self.analyzer._score_mhc_ii_beta(beta_seq, len(beta_seq))
        self.assertGreater(beta_score, 0.5)
    
    def test_tcr_scoring(self):
        """Test TCR scoring functions (should be same as MHC-I)."""
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
        self.assertIn('length_range', self.analyzer.mhc_ii_patterns['alpha'])
        self.assertIn('length_range', self.analyzer.mhc_ii_patterns['beta'])
        self.assertIn('highly_specific_motifs', self.analyzer.mhc_ii_patterns['alpha'])
        self.assertIn('highly_specific_motifs', self.analyzer.mhc_ii_patterns['beta'])
        
        # Check MHC-II specific patterns
        alpha_motifs = self.analyzer.mhc_ii_patterns['alpha']['highly_specific_motifs']
        self.assertIn('MGSGWV', alpha_motifs)  # MHC-II alpha signature
        
        beta_motifs = self.analyzer.mhc_ii_patterns['beta']['highly_specific_motifs']
        self.assertIn('MSWKKA', beta_motifs)  # MHC-II beta signature
    
    def test_file_not_found(self):
        """Test handling of non-existent files."""
        with self.assertRaises(FileNotFoundError):
            self.analyzer.analyze_pdb("nonexistent_file.pdb")
    
    def create_mock_mhc_ii_pdb(self, complex_type="simple") -> str:
        """Create mock PDB files for testing MHC-II scenarios."""
        
        if complex_type == "simple":
            # Simple MHC-II complex with 2 chains
            pdb_content = """HEADER    IMMUNE SYSTEM/IMMUNE SYSTEM             01-JAN-20   TEST    
ATOM      1  N   MET A   1      20.154  16.000  12.000  1.00 10.00           N  
ATOM      2  CA  MET A   1      19.000  16.500  12.500  1.00 10.00           C  
ATOM      3  N   GLY B   1      25.000  20.000  15.000  1.00 10.00           N  
ATOM      4  CA  GLY B   1      24.500  20.500  15.500  1.00 10.00           C  
END
"""
        elif complex_type == "mhc_ii_complex":
            # Mock MHC-II complex with realistic chain lengths
            pdb_content = """HEADER    IMMUNE SYSTEM/IMMUNE SYSTEM             01-JAN-20   MOCK    
TITLE     MOCK MHC-II-TCR COMPLEX
"""
            # Add longer peptide chain (chain P) - 16 residues (MHC-II length)
            atom_num = 1
            peptide_aas = ['TYR', 'LEU', 'GLN', 'PRO', 'ARG', 'TRP', 'TYR', 'PHE', 'VAL', 
                          'ALA', 'SER', 'THR', 'GLU', 'ASP', 'LYS', 'ARG']
            for i, aa in enumerate(peptide_aas):
                x, y, z = 20.0 + i, 15.0, 10.0
                pdb_content += f"ATOM  {atom_num:5d}  N   {aa} P{i+1:4d}    {x:8.3f}{y:8.3f}{z:8.3f}  1.00 20.00           N  \n"
                atom_num += 1
                pdb_content += f"ATOM  {atom_num:5d}  CA  {aa} P{i+1:4d}    {x+0.5:8.3f}{y+0.5:8.3f}{z+0.5:8.3f}  1.00 20.00           C  \n"
                atom_num += 1
            
            # Add MHC-II alpha chain (chain A) - 190 residues (simplified)
            for i in range(190):
                aa = ['ALA', 'VAL', 'LEU', 'ILE'][i % 4]
                x, y, z = 30.0 + i*0.1, 20.0, 15.0
                pdb_content += f"ATOM  {atom_num:5d}  N   {aa} A{i+1:4d}    {x:8.3f}{y:8.3f}{z:8.3f}  1.00 20.00           N  \n"
                atom_num += 1
                pdb_content += f"ATOM  {atom_num:5d}  CA  {aa} A{i+1:4d}    {x+0.5:8.3f}{y+0.5:8.3f}{z+0.5:8.3f}  1.00 20.00           C  \n"
                atom_num += 1
            
            # Add MHC-II beta chain (chain B) - 200 residues
            for i in range(200):
                aa = ['SER', 'THR', 'ASN', 'GLN'][i % 4]
                x, y, z = 40.0 + i*0.1, 25.0, 20.0
                pdb_content += f"ATOM  {atom_num:5d}  N   {aa} B{i+1:4d}    {x:8.3f}{y:8.3f}{z:8.3f}  1.00 20.00           N  \n"
                atom_num += 1
                pdb_content += f"ATOM  {atom_num:5d}  CA  {aa} B{i+1:4d}    {x+0.5:8.3f}{y+0.5:8.3f}{z+0.5:8.3f}  1.00 20.00           C  \n"
                atom_num += 1
            
            # Add TCR alpha chain (chain T) - 210 residues
            for i in range(210):
                aa = ['GLU', 'ASP', 'LYS', 'ARG'][i % 4]
                x, y, z = 50.0 + i*0.1, 30.0, 25.0
                pdb_content += f"ATOM  {atom_num:5d}  N   {aa} T{i+1:4d}    {x:8.3f}{y:8.3f}{z:8.3f}  1.00 20.00           N  \n"
                atom_num += 1
                pdb_content += f"ATOM  {atom_num:5d}  CA  {aa} T{i+1:4d}    {x+0.5:8.3f}{y+0.5:8.3f}{z+0.5:8.3f}  1.00 20.00           C  \n"
                atom_num += 1
                
            # Add TCR beta chain (chain U) - 250 residues
            for i in range(250):
                aa = ['PHE', 'TYR', 'TRP', 'HIS'][i % 4]
                x, y, z = 60.0 + i*0.1, 35.0, 30.0
                pdb_content += f"ATOM  {atom_num:5d}  N   {aa} U{i+1:4d}    {x:8.3f}{y:8.3f}{z:8.3f}  1.00 20.00           N  \n"
                atom_num += 1
                pdb_content += f"ATOM  {atom_num:5d}  CA  {aa} U{i+1:4d}    {x+0.5:8.3f}{y+0.5:8.3f}{z+0.5:8.3f}  1.00 20.00           C  \n"
                atom_num += 1
            
            pdb_content += "END\n"
        
        # Create temporary file
        temp_file = tempfile.NamedTemporaryFile(mode='w', suffix='.pdb', delete=False)
        temp_file.write(pdb_content)
        temp_file.close()
        
        return temp_file.name
    
    def test_pdb_parsing(self):
        """Test basic PDB parsing functionality."""
        mock_pdb = self.create_mock_mhc_ii_pdb("simple")
        
        try:
            result = self.analyzer.analyze_pdb(mock_pdb)
            self.assertIsInstance(result, dict)
            self.assertGreater(len(result), 0)
        finally:
            os.unlink(mock_pdb)
    
    def test_mhc_ii_complex(self):
        """Test analysis of a mock MHC-II-TCR complex."""
        mock_pdb = self.create_mock_mhc_ii_pdb("mhc_ii_complex")
        
        try:
            result = self.analyzer.analyze_pdb(mock_pdb)
            self.assertIsInstance(result, dict)
            
            # Should have 5 chains
            self.assertEqual(len(result), 5)
            
            # Check that we have all expected chain types (with or without complex numbers)
            chain_types = set()
            for chain_type in result.values():
                # Remove complex numbers (e.g., "peptide_complex1" -> "peptide")
                base_type = chain_type.split('_complex')[0] if '_complex' in chain_type else chain_type
                chain_types.add(base_type)
            
            expected_types = {'peptide', 'mhc_ii_alpha', 'mhc_ii_beta', 'tcr_alpha', 'tcr_beta'}
            
            # At least 3 of the 5 expected types should be identified
            overlap = len(chain_types.intersection(expected_types))
            self.assertGreaterEqual(overlap, 3)
            
            if self.analyzer.verbose:
                print(f"Mock MHC-II complex analysis result: {result}")
            
        finally:
            os.unlink(mock_pdb)
    
    def test_chain_classification(self):
        """Test chain classification logic."""
        # Create mock chain info
        chain_info = {
            'P': {'sequence': 'PKYVKQNTLKLAT', 'length': 13, 'properties': {'molecular_weight': 1500}},
            'A': {'sequence': 'MGSGWV' + 'A' * 180, 'length': 186, 'properties': {'molecular_weight': 20000}},
            'B': {'sequence': 'MSWKKA' + 'B' * 190, 'length': 196, 'properties': {'molecular_weight': 21000}},
        }
        
        assignments = self.analyzer._classify_chains(chain_info)
        
        # Check that assignment was made
        self.assertEqual(len(assignments), 3)
        self.assertIn('P', assignments)
        self.assertIn('A', assignments)
        self.assertIn('B', assignments)
    
    def test_scoring_functions(self):
        """Test individual scoring functions."""
        # Test peptide scoring
        peptide_score = self.analyzer._score_peptide("PKYVKQNTLKLAT", 13, {'molecular_weight': 1500})
        self.assertGreater(peptide_score, 0.0)
        
        # Test MHC-II alpha scoring
        alpha_score = self.analyzer._score_mhc_ii_alpha("MGSGWV" + "A" * 180, 186)
        self.assertGreater(alpha_score, 0.0)
        
        # Test MHC-II beta scoring
        beta_score = self.analyzer._score_mhc_ii_beta("MSWKKA" + "B" * 190, 196)
        self.assertGreater(beta_score, 0.0)


if __name__ == '__main__':
    unittest.main()