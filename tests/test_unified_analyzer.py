#!/usr/bin/env python3
"""
Unit tests for the Unified TCR-pMHC Analyzer

Tests both MHC-I and MHC-II complex analysis, chain identification,
clustering, and pairing functionality.
"""

import unittest
import tempfile
import os
from pathlib import Path
import numpy as np

# Import the unified analyzer
from phlatcr_splicer import TCRpMHCAnalyzer


class TestTCRpMHCAnalyzer(unittest.TestCase):
    """Test cases for the unified TCR-pMHC analyzer."""
    
    def setUp(self):
        """Set up test fixtures."""
        self.analyzer = TCRpMHCAnalyzer(verbose=False)
        self.test_data_dir = Path(__file__).parent / "test_data"
    
    def test_initialization(self):
        """Test analyzer initialization with different parameters."""
        # Default initialization
        analyzer = TCRpMHCAnalyzer()
        self.assertEqual(analyzer.eps, 60.0)
        self.assertEqual(analyzer.align_score, 20)
        
        # Custom parameters
        analyzer = TCRpMHCAnalyzer(
            eps=45.0,
            align_score=25,
            tcr_pair_dist=40.0
        )
        self.assertEqual(analyzer.eps, 45.0)
        self.assertEqual(analyzer.align_score, 25)
        self.assertEqual(analyzer.tcr_pair_dist, 40.0)
    
    def test_chain_sequence_extraction(self):
        """Test extraction of amino acid sequences from chains."""
        # Create a mock chain structure
        sequence = "ACDEFGHIKLMNPQRSTVWY"
        extracted = self.analyzer._get_chain_sequence(self._create_mock_chain(sequence))
        self.assertEqual(extracted, sequence)
    
    def test_peptide_identification(self):
        """Test peptide identification by length."""
        # Short peptide (MHC-I range)
        short_peptide = "NLVPMVATV"  # 9 residues
        chain_type = self.analyzer._identify_chain_type(short_peptide)
        self.assertEqual(chain_type, 'PEPTIDE')
        
        # Long peptide (MHC-II range)
        long_peptide = "PKYVKQNTLKLATGMR"  # 16 residues
        chain_type = self.analyzer._identify_chain_type(long_peptide)
        self.assertEqual(chain_type, 'PEPTIDE')
        
        # Edge case: 40 residues (max peptide length)
        max_peptide = "A" * 40
        chain_type = self.analyzer._identify_chain_type(max_peptide)
        self.assertEqual(chain_type, 'PEPTIDE')
        
        # Too long for peptide
        not_peptide = "A" * 41
        chain_type = self.analyzer._identify_chain_type(not_peptide)
        self.assertNotEqual(chain_type, 'PEPTIDE')
    
    def test_tcr_identification_by_alignment(self):
        """Test TCR alpha/beta discrimination using sequence alignment."""
        # Create analyzer with lower alignment threshold for testing
        analyzer = TCRpMHCAnalyzer(align_score=10)
        
        # Test with a sequence containing TCR-alpha motifs
        tcr_alpha_like = ("MAMLELLTLAFLGIWAFNQRADGQIYNQEPATNENISEATYNASLCSLTESGKSYFFWYKQEP"
                         "GAGLQLLTYIFSNMDMKQDQRLTVLLLNKKDKHLSLRIADTQTGDSAIYFCAVSGGYQKVTFG"
                         "IGTKLQVIP")
        chain_type = analyzer._identify_chain_type(tcr_alpha_like)
        # Should identify as TCR based on length and patterns
        self.assertIn(chain_type, ['TCR_ALPHA', 'TCR_BETA'])
    
    def test_mhc_identification_by_pattern(self):
        """Test MHC chain identification using patterns."""
        # MHC-I heavy chain pattern
        mhc_i_seq = "GSHSMR" + "A" * 260  # Contains MHC-I specific pattern
        chain_type = self.analyzer._identify_chain_type(mhc_i_seq)
        self.assertEqual(chain_type, 'MHC_I_ALPHA')
        
        # B2M pattern
        b2m_seq = "IQRTPKIQ" + "A" * 90  # Contains B2M specific pattern
        chain_type = self.analyzer._identify_chain_type(b2m_seq)
        self.assertEqual(chain_type, 'B2M')
        
        # MHC-II alpha pattern
        mhc_ii_alpha_seq = "WRLEEF" + "A" * 180  # Contains MHC-II alpha pattern
        chain_type = self.analyzer._identify_chain_type(mhc_ii_alpha_seq)
        self.assertEqual(chain_type, 'MHC_II_ALPHA')
        
        # MHC-II beta pattern
        mhc_ii_beta_seq = "NGTER" + "A" * 190  # Contains MHC-II beta pattern
        chain_type = self.analyzer._identify_chain_type(mhc_ii_beta_seq)
        self.assertEqual(chain_type, 'MHC_II_BETA')
    
    def test_chain_center_computation(self):
        """Test computation of chain geometric centers."""
        # Create a mock chain with known coordinates
        mock_chain = self._create_mock_chain_with_coords([
            [0, 0, 0],
            [1, 1, 1],
            [2, 2, 2]
        ])
        center = self.analyzer._compute_chain_center(mock_chain)
        expected_center = np.array([1, 1, 1])  # Mean of the coordinates
        np.testing.assert_array_almost_equal(center, expected_center)
    
    def test_unknown_chain_reclassification(self):
        """Test reclassification of UNKNOWN chains based on proximity."""
        # Set up chain types and centers
        chain_types = {
            'A': 'MHC_I_ALPHA',
            'B': 'UNKNOWN',  # Should become B2M
            'C': 'MHC_II_BETA',
            'D': 'UNKNOWN',  # Should become MHC_II_ALPHA
        }
        
        chain_centers = {
            'A': np.array([0, 0, 0]),
            'B': np.array([10, 0, 0]),  # Close to MHC_I_ALPHA
            'C': np.array([100, 0, 0]),
            'D': np.array([110, 0, 0]),  # Close to MHC_II_BETA
        }
        
        # Create mock structure
        mock_structure = self._create_mock_structure_for_reclassification()
        
        # Run reclassification
        self.analyzer._reclassify_unknown_chains(chain_types, chain_centers, mock_structure)
        
        # Check reclassification
        self.assertEqual(chain_types['B'], 'B2M')
        self.assertEqual(chain_types['D'], 'MHC_II_ALPHA')
    
    def test_complex_pairing_mhc_i(self):
        """Test pairing of chains in MHC-I complexes."""
        chains = ['A', 'B', 'C', 'D', 'E']
        types = {
            'A': 'MHC_I_ALPHA',
            'B': 'B2M',
            'C': 'PEPTIDE',
            'D': 'TCR_ALPHA',
            'E': 'TCR_BETA'
        }
        centers = {
            'A': np.array([0, 0, 0]),
            'B': np.array([5, 0, 0]),    # Close to MHC_I_ALPHA
            'C': np.array([2.5, 0, 0]),  # Between MHC_I_ALPHA and B2M
            'D': np.array([20, 0, 0]),   # TCR alpha
            'E': np.array([25, 0, 0]),   # TCR beta, close to alpha
        }
        
        mock_structure = self._create_mock_structure_for_reclassification()
        pairs = self.analyzer._pair_chains_in_group(mock_structure, chains, types, centers)
        
        # Should form one complete TCR-pMHC-I complex
        self.assertEqual(len(pairs), 1)
        complex_pair = pairs[0]
        self.assertEqual(complex_pair['type'], 'MHC-I')
        self.assertEqual(complex_pair['mhc_alpha'], 'A')
        self.assertEqual(complex_pair['b2m'], 'B')
        self.assertEqual(complex_pair['peptide'], 'C')
        self.assertEqual(complex_pair['tcr_alpha'], 'D')
        self.assertEqual(complex_pair['tcr_beta'], 'E')
    
    def test_complex_pairing_mhc_ii(self):
        """Test pairing of chains in MHC-II complexes."""
        chains = ['A', 'B', 'C', 'D', 'E']
        types = {
            'A': 'MHC_II_ALPHA',
            'B': 'MHC_II_BETA',
            'C': 'PEPTIDE',
            'D': 'TCR_ALPHA',
            'E': 'TCR_BETA'
        }
        centers = {
            'A': np.array([0, 0, 0]),
            'B': np.array([5, 0, 0]),    # Close to MHC_II_ALPHA
            'C': np.array([2.5, 0, 0]),  # Between MHC_II chains
            'D': np.array([20, 0, 0]),   # TCR alpha
            'E': np.array([25, 0, 0]),   # TCR beta, close to alpha
        }
        
        mock_structure = self._create_mock_structure_for_reclassification()
        pairs = self.analyzer._pair_chains_in_group(mock_structure, chains, types, centers)
        
        # Should form one complete TCR-pMHC-II complex
        self.assertEqual(len(pairs), 1)
        complex_pair = pairs[0]
        self.assertEqual(complex_pair['type'], 'MHC-II')
        self.assertEqual(complex_pair['mhc_alpha'], 'A')
        self.assertEqual(complex_pair['mhc_beta'], 'B')
        self.assertEqual(complex_pair['peptide'], 'C')
        self.assertEqual(complex_pair['tcr_alpha'], 'D')
        self.assertEqual(complex_pair['tcr_beta'], 'E')
    
    def test_unpaired_tcr_handling(self):
        """Test handling of unpaired TCR chains."""
        chains = ['A', 'B']
        types = {
            'A': 'TCR_ALPHA',
            'B': 'TCR_BETA'
        }
        centers = {
            'A': np.array([0, 0, 0]),
            'B': np.array([5, 0, 0]),  # Close to TCR_ALPHA
        }
        
        mock_structure = self._create_mock_structure_for_reclassification()
        pairs = self.analyzer._pair_chains_in_group(mock_structure, chains, types, centers)
        
        # Should form one unpaired TCR complex
        self.assertEqual(len(pairs), 1)
        complex_pair = pairs[0]
        self.assertEqual(complex_pair['type'], 'UNPAIRED_TCR')
        self.assertEqual(complex_pair['tcr_alpha'], 'A')
        self.assertEqual(complex_pair['tcr_beta'], 'B')
    
    def test_file_format_support(self):
        """Test that both PDB and CIF files are supported."""
        # Test with a PDB file path
        pdb_file = "test.pdb"
        # Would need actual file to test fully
        
        # Test with a CIF file path
        cif_file = "test.cif"
        # Would need actual file to test fully
    
    # Helper methods for creating mock structures
    
    def _create_mock_chain(self, sequence):
        """Create a mock chain with the given sequence."""
        class MockResidue:
            def __init__(self, resname):
                self.resname = resname
                self.id = (' ', 1, ' ')
        
        class MockChain:
            def __init__(self, sequence):
                # Convert single letter to three letter codes
                aa_map = {
                    'A': 'ALA', 'C': 'CYS', 'D': 'ASP', 'E': 'GLU',
                    'F': 'PHE', 'G': 'GLY', 'H': 'HIS', 'I': 'ILE',
                    'K': 'LYS', 'L': 'LEU', 'M': 'MET', 'N': 'ASN',
                    'P': 'PRO', 'Q': 'GLN', 'R': 'ARG', 'S': 'SER',
                    'T': 'THR', 'V': 'VAL', 'W': 'TRP', 'Y': 'TYR'
                }
                self.residues = [MockResidue(aa_map.get(aa, 'UNK')) for aa in sequence]
            
            def __iter__(self):
                return iter(self.residues)
        
        return MockChain(sequence)
    
    def _create_mock_chain_with_coords(self, coords):
        """Create a mock chain with specific CA atom coordinates."""
        class MockAtom:
            def __init__(self, coord):
                self.name = 'CA'
                self.coord = np.array(coord)
        
        class MockResidue:
            def __init__(self, coord):
                self.id = (' ', 1, ' ')
                self.atoms = [MockAtom(coord)]
            
            def __iter__(self):
                return iter(self.atoms)
        
        class MockChain:
            def __init__(self, coords):
                self.residues = [MockResidue(coord) for coord in coords]
            
            def __iter__(self):
                return iter(self.residues)
        
        return MockChain(coords)
    
    def _create_mock_structure_for_reclassification(self):
        """Create a mock structure for reclassification tests."""
        class MockModel:
            def __init__(self):
                self.chains = {
                    'A': self._create_mock_chain('A' * 270),  # MHC_I_ALPHA length
                    'B': self._create_mock_chain('A' * 100),  # B2M length
                    'C': self._create_mock_chain('A' * 200),  # MHC_II_BETA length
                    'D': self._create_mock_chain('A' * 180),  # MHC_II_ALPHA length
                }
            
            def __getitem__(self, key):
                return self.chains[key]
        
        class MockStructure:
            def __init__(self):
                self.models = [MockModel()]
            
            def __getitem__(self, idx):
                return self.models[idx]
        
        # Need to reference the test class's helper method
        test_instance = TestTCRpMHCAnalyzer()
        MockModel._create_mock_chain = test_instance._create_mock_chain
        
        return MockStructure()


class TestIntegration(unittest.TestCase):
    """Integration tests for the unified analyzer."""
    
    def setUp(self):
        """Set up test fixtures."""
        self.analyzer = TCRpMHCAnalyzer(verbose=False)
    
    def test_create_temp_pdb_file(self):
        """Test analysis of a temporary PDB file."""
        # Create a minimal PDB file
        pdb_content = """HEADER    TEST STRUCTURE
ATOM      1  CA  ALA A   1       0.000   0.000   0.000  1.00  0.00           C
ATOM      2  CA  ALA A   2       3.800   0.000   0.000  1.00  0.00           C
ATOM      3  CA  ALA A   3       7.600   0.000   0.000  1.00  0.00           C
END
"""
        with tempfile.NamedTemporaryFile(mode='w', suffix='.pdb', delete=False) as f:
            f.write(pdb_content)
            temp_file = f.name
        
        try:
            # Analyze the file
            results, chain_map = self.analyzer.analyze_pdb(temp_file)
            # Should process without errors even if no valid complexes found
            self.assertIsInstance(results, dict)
            self.assertIsInstance(chain_map, dict)
        finally:
            # Clean up
            os.unlink(temp_file)


if __name__ == '__main__':
    unittest.main()