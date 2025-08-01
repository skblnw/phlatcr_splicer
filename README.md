# pHLA-TCR Complex Structure Analyzer

A Python tool for analyzing protein complex structure PDB files, specifically designed for pHLA-TCR complexes. This tool identifies and classifies protein chains even when chain IDs are inconsistent.

## Features

- **Chain Identification**: Automatically identifies chains in pHLA-TCR complexes:
  - Peptide antigen
  - MHC heavy chain
  - β2-microglobulin (b2m)
  - TCR alpha chain
  - TCR beta chain

- **Robust Analysis**: Works with inconsistent chain IDs by analyzing structural and sequence features

- **PDB Support**: Handles standard PDB file formats

## Installation

```bash
git clone https://github.com/yourusername/phlatcr_splicer.git
cd phlatcr_splicer
pip install -r requirements.txt
```

## Usage

```python
from phlatcr_analyzer import pHLATCRAnalyzer

# Initialize analyzer
analyzer = pHLATCRAnalyzer()

# Analyze a PDB file
result = analyzer.analyze_pdb("complex.pdb")

# Print chain assignments
for chain_id, chain_type in result.items():
    print(f"Chain {chain_id}: {chain_type}")
```

## Requirements

- Python 3.7+
- BioPython
- NumPy
- SciPy

## License

BSD 3-Clause License

## Citation

If you use this tool in your research, please cite:
[Citation information to be added]