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

if __name__ == '__main__':
    unittest.main()