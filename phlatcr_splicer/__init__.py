"""
pHLA-TCR Complex Structure Analyzer

A Python package for analyzing protein complex structure PDB files, 
specifically designed for pHLA-TCR complexes. This tool identifies 
and classifies protein chains even when chain IDs are inconsistent.
"""

from .analyzer import pHLATCRAnalyzer

__version__ = "0.1.0"
__author__ = "Your Name"
__email__ = "your.email@example.com"

__all__ = ["pHLATCRAnalyzer"]