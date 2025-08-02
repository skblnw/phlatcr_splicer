"""
pHLA-TCR Complex Structure Analyzer

A Python package for analyzing protein complex structure PDB files, 
specifically designed for pHLA-TCR complexes. This tool identifies 
and classifies protein chains even when chain IDs are inconsistent.
"""

from .analyzer import pMHCITCRAnalyzer
from .mhc_ii_analyzer import pMHCIITCRAnalyzer

__all__ = ['pMHCITCRAnalyzer', 'pMHCIITCRAnalyzer']

__version__ = "0.1.0"
__author__ = "Your Name"
__email__ = "your.email@example.com"

__all__ = ["pHLATCRAnalyzer"]