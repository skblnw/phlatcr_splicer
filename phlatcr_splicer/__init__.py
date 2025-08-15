"""
pHLA-TCR Complex Structure Analyzer

A unified Python package for analyzing TCR-pMHC protein complex structures
from PDB and CIF files. Automatically handles both MHC-I and MHC-II complexes
using sequence alignment, spatial clustering, and pattern recognition.
"""

# Import unified analyzer only
from .unified_analyzer import TCRpMHCAnalyzer

__all__ = ['TCRpMHCAnalyzer']

__version__ = "1.0.0"
__author__ = "skblnw"
__email__ = "skblnw@github.com"