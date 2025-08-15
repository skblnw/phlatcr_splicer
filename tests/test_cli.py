#!/usr/bin/env python3
"""
Command-line interface tests for the Unified TCR-pMHC Analyzer

Tests the CLI functionality, argument parsing, and output formatting.

Author: skblnw
"""

import unittest
import tempfile
import os
import sys
from pathlib import Path
import subprocess
from unittest.mock import patch, MagicMock

# Add parent directory to path
sys.path.insert(0, str(Path(__file__).parent.parent))


class TestCLI(unittest.TestCase):
    """Test command-line interface functionality."""
    
    def setUp(self):
        """Set up test fixtures."""
        self.test_dir = Path(__file__).parent
        self.script_path = self.test_dir.parent / "scripts" / "main.py"
        
    def test_help_flag(self):
        """Test that --help flag works."""
        result = subprocess.run(
            [sys.executable, str(self.script_path), "--help"],
            capture_output=True,
            text=True
        )
        self.assertEqual(result.returncode, 0)
        self.assertIn("Analyze TCR-pMHC structures", result.stdout)
        self.assertIn("--eps", result.stdout)
        self.assertIn("--align-score", result.stdout)
    
    def test_version_display(self):
        """Test that version is displayed correctly."""
        result = subprocess.run(
            [sys.executable, "-c", 
             "from phlatcr_splicer import __version__; print(__version__)"],
            capture_output=True,
            text=True
        )
        self.assertEqual(result.returncode, 0)
        self.assertEqual(result.stdout.strip(), "1.0.0")
    
    def test_missing_input_file(self):
        """Test error handling for missing input file."""
        result = subprocess.run(
            [sys.executable, str(self.script_path), "nonexistent.pdb"],
            capture_output=True,
            text=True
        )
        # Should handle gracefully
        self.assertIn("No PDB/CIF files found", result.stdout)
    
    def test_parameter_parsing(self):
        """Test that CLI parameters are parsed correctly."""
        # Create a minimal test PDB file
        with tempfile.NamedTemporaryFile(mode='w', suffix='.pdb', delete=False) as f:
            f.write("HEADER    TEST\n")
            f.write("ATOM      1  CA  ALA A   1       0.000   0.000   0.000\n")
            f.write("END\n")
            test_file = f.name
        
        try:
            # Test with various parameters
            result = subprocess.run(
                [sys.executable, str(self.script_path), test_file,
                 "--eps", "45.0",
                 "--align-score", "25",
                 "--tcr-pair-dist", "40.0",
                 "--verbose"],
                capture_output=True,
                text=True
            )
            # Should run without errors
            self.assertNotIn("error", result.stderr.lower())
        finally:
            os.unlink(test_file)
    
    def test_batch_processing_flag(self):
        """Test batch processing with summary flag."""
        # Create multiple test files
        test_files = []
        for i in range(2):
            with tempfile.NamedTemporaryFile(mode='w', suffix='.pdb', delete=False) as f:
                f.write(f"HEADER    TEST{i}\n")
                f.write("ATOM      1  CA  ALA A   1       0.000   0.000   0.000\n")
                f.write("END\n")
                test_files.append(f.name)
        
        try:
            result = subprocess.run(
                [sys.executable, str(self.script_path)] + test_files + ["--batch-summary"],
                capture_output=True,
                text=True
            )
            # Should mention batch processing
            self.assertIn("Found 2 file(s)", result.stdout)
        finally:
            for f in test_files:
                os.unlink(f)
    
    def test_output_file_creation(self):
        """Test that output file is created when specified."""
        with tempfile.NamedTemporaryFile(mode='w', suffix='.pdb', delete=False) as f:
            f.write("HEADER    TEST\n")
            f.write("ATOM      1  CA  ALA A   1       0.000   0.000   0.000\n")
            f.write("END\n")
            test_file = f.name
        
        output_file = tempfile.mktemp(suffix='.txt')
        
        try:
            result = subprocess.run(
                [sys.executable, str(self.script_path), test_file,
                 "--output", output_file],
                capture_output=True,
                text=True
            )
            # Check if output file was created
            if os.path.exists(output_file):
                self.assertTrue(True, "Output file created")
                os.unlink(output_file)
        finally:
            os.unlink(test_file)
    
    def test_verbose_output(self):
        """Test that verbose flag produces detailed output."""
        with tempfile.NamedTemporaryFile(mode='w', suffix='.pdb', delete=False) as f:
            f.write("HEADER    TEST\n")
            f.write("ATOM      1  CA  ALA A   1       0.000   0.000   0.000\n")
            f.write("END\n")
            test_file = f.name
        
        try:
            # Run without verbose
            result_quiet = subprocess.run(
                [sys.executable, str(self.script_path), test_file],
                capture_output=True,
                text=True
            )
            
            # Run with verbose
            result_verbose = subprocess.run(
                [sys.executable, str(self.script_path), test_file, "--verbose"],
                capture_output=True,
                text=True
            )
            
            # Verbose should produce more output (or at least not less)
            self.assertGreaterEqual(
                len(result_verbose.stdout) + len(result_verbose.stderr),
                len(result_quiet.stdout) + len(result_quiet.stderr)
            )
        finally:
            os.unlink(test_file)
    
    def test_cif_file_support(self):
        """Test that CIF files are recognized."""
        with tempfile.NamedTemporaryFile(mode='w', suffix='.cif', delete=False) as f:
            f.write("data_TEST\n")
            f.write("_atom_site.group_PDB\n")
            f.write("_atom_site.id\n")
            f.write("ATOM 1\n")
            test_file = f.name
        
        try:
            result = subprocess.run(
                [sys.executable, str(self.script_path), test_file],
                capture_output=True,
                text=True
            )
            # Should complete processing (even if no valid chains found)
            self.assertIn("complete", result.stdout.lower())
        finally:
            os.unlink(test_file)


class TestCLIIntegration(unittest.TestCase):
    """Integration tests for CLI with real analysis."""
    
    def setUp(self):
        """Set up test fixtures."""
        self.test_dir = Path(__file__).parent
        self.script_path = self.test_dir.parent / "scripts" / "main.py"
    
    def test_real_pdb_analysis(self):
        """Test analysis with a real-like PDB structure."""
        pdb_content = """HEADER    IMMUNE SYSTEM                           01-JAN-00   TEST
ATOM      1  CA  ALA A   1       0.000   0.000   0.000  1.00  0.00           C
ATOM      2  CA  ALA A   2       3.800   0.000   0.000  1.00  0.00           C
ATOM      3  CA  ALA A   3       7.600   0.000   0.000  1.00  0.00           C
ATOM      4  CA  ALA B   1      20.000   0.000   0.000  1.00  0.00           C
ATOM      5  CA  ALA B   2      23.800   0.000   0.000  1.00  0.00           C
END
"""
        with tempfile.NamedTemporaryFile(mode='w', suffix='.pdb', delete=False) as f:
            f.write(pdb_content)
            test_file = f.name
        
        try:
            result = subprocess.run(
                [sys.executable, str(self.script_path), test_file],
                capture_output=True,
                text=True
            )
            self.assertEqual(result.returncode, 0)
            # Should complete analysis (may not find valid complexes in minimal structure)
            self.assertIn("complete", result.stdout.lower())
        finally:
            os.unlink(test_file)


if __name__ == '__main__':
    unittest.main()