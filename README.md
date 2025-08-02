# pHLA-TCR Splicer

A comprehensive Python tool for analyzing protein complex structures from PDB files, supporting both **MHC Class I** and **MHC Class II** complexes. The tool automatically identifies and classifies different chain types even when chain IDs are inconsistent across PDB files.

## Features

- **Dual Platform Support**: Analyze both MHC-I (pHLA-TCR) and MHC-II (pMHC-II-TCR) complexes
- **Smart Auto-Detection**: Automatically detect complex type (MHC-I vs MHC-II)
- **Batch Processing**: Analyze multiple PDB files at once with summary statistics
- **Advanced Pattern Recognition**: Uses sequence motifs, length analysis, and molecular properties
- **Multi-Complex Support**: Handles PDB files with multiple complexes
- **Flexible Interface**: Command-line, Python API, and individual analyzer access
- **Detailed Reporting**: Confidence scores and comprehensive analysis reports

## Quick Start

### Installation

```bash
git clone https://github.com/yourusername/phlatcr_splicer.git
cd phlatcr_splicer
pip install -e .
```

### Basic Usage

```bash
# Single file analysis
python scripts/main.py --type mhc-i test_data/1oga.pdb
python scripts/main.py --type mhc-ii test_data/4z7u.pdb
python scripts/main.py --auto test_data/complex.pdb --verbose

# Batch processing (multiple files)
python scripts/main.py --type mhc-i *.pdb --batch-summary
python scripts/main.py --auto file1.pdb file2.pdb file3.pdb --output results.txt

# Individual analyzers
phlatcr-analyze test_data/1oga.pdb          # MHC-I specific
mhc-ii-analyze test_data/4z7u.pdb           # MHC-II specific
```

## Usage Options

### 1. Unified Main Script (Recommended)

```bash
# Single file analysis
python scripts/main.py --type mhc-i input.pdb
python scripts/main.py --type mhc-ii input.pdb
python scripts/main.py --auto input.pdb --verbose

# Batch processing multiple files
python scripts/main.py --type mhc-i file1.pdb file2.pdb file3.pdb
python scripts/main.py --auto *.pdb --batch-summary
python scripts/main.py --type mhc-ii complexes/*.pdb --output batch_results.txt

# Save results to file
python scripts/main.py --type mhc-i input.pdb --output results.txt

# Help
python scripts/main.py --help
```

### 2. Individual Console Commands

```bash
# MHC-I analyzer
phlatcr-analyze input.pdb --verbose

# MHC-II analyzer  
mhc-ii-analyze input.pdb --verbose
```

### 3. Python API

```python
from phlatcr_splicer import pMHCITCRAnalyzer, pMHCIITCRAnalyzer

# MHC-I Analysis
mhc_i_analyzer = pMHCITCRAnalyzer(verbose=True)
mhc_i_results = mhc_i_analyzer.analyze_pdb("test_data/1oga.pdb")

# MHC-II Analysis
mhc_ii_analyzer = pMHCIITCRAnalyzer(verbose=True)
mhc_ii_results = mhc_ii_analyzer.analyze_pdb("test_data/4z7u.pdb")

print(mhc_i_results)
# {'A': 'mhc_heavy', 'B': 'b2m', 'C': 'peptide', 'D': 'tcr_alpha', 'E': 'tcr_beta'}
```

## Supported Complex Types

### MHC Class I (pHLA-TCR)
- **`mhc_heavy`**: MHC heavy chain (~270-380 residues)
- **`b2m`**: β2-microglobulin (~99 residues)
- **`peptide`**: Short peptides (8-11 residues)
- **`tcr_alpha`**: TCR alpha chain (~180-230 residues)  
- **`tcr_beta`**: TCR beta chain (~230-290 residues)

### MHC Class II (pMHC-II-TCR)
- **`mhc_ii_alpha`**: MHC-II alpha chain (~180-200 residues)
- **`mhc_ii_beta`**: MHC-II beta chain (~190-210 residues)
- **`peptide`**: Longer peptides (12-25 residues)
- **`tcr_alpha`**: TCR alpha chain (~180-230 residues)
- **`tcr_beta`**: TCR beta chain (~230-290 residues)

## Example Results

### MHC-I Complex (1oga.pdb)
```bash
python scripts/main.py --type mhc-i test_data/1oga.pdb
```
```
Running MHC-I (pHLA-TCR) Analysis
==================================================

Analysis Results for 1oga.pdb:
----------------------------------------
  Chain A: mhc_heavy
  Chain B: b2m
  Chain C: peptide
  Chain D: tcr_alpha
  Chain E: tcr_beta

Summary:
  Total chains: 5
  Identified: 5
  Unknown: 0
Perfect! All chains identified successfully.
```

### MHC-II Multi-Complex (4z7u.pdb)
```bash
python scripts/main.py --type mhc-ii test_data/4z7u.pdb
```
```
Running MHC-II (pMHC-II-TCR) Analysis
==================================================

Analysis Results for 4z7u.pdb:
----------------------------------------
  Chain A: mhc_ii_alpha_complex1
  Chain B: mhc_ii_beta_complex1
  Chain C: mhc_ii_alpha_complex2
  Chain D: mhc_ii_beta_complex2
  Chain E: tcr_alpha_complex1
  Chain F: tcr_beta_complex1
  Chain G: tcr_alpha_complex2
  Chain H: tcr_beta_complex2
  Chain I: peptide_complex1
  Chain J: peptide_complex2

Summary:
  Total chains: 10
  Identified: 10
  Unknown: 0
Perfect! All chains identified successfully.
```

### Batch Processing

```bash
python scripts/main.py --type mhc-i test_data/1oga.pdb test_data/1mi5.pdb --batch-summary
```
```
Running Batch Analysis
==================================================
Files to process: 2

Processing 1/2: 1oga.pdb
   ✅ 5/5 chains identified

Processing 2/2: 1mi5.pdb  
   ✅ 5/5 chains identified

Batch Analysis Summary
==============================
Total files processed: 2
Successful analyses: 2
Failed analyses: 0
Total chains: 10
Identified chains: 10
Unknown chains: 0
Identification rate: 100.0%

Batch analysis completed successfully!
```

### Auto-Detection
```bash
python scripts/main.py --auto test_data/1oga.pdb --verbose
```
```
Auto-detecting complex type...
Detected: MHC-I complex (β2-microglobulin detected)
   MHC-I analysis: 5/5 chains identified
   MHC-II analysis: 2/5 chains identified
Running MHC-I (pHLA-TCR) Analysis
...
```

## Testing

```bash
# Run all tests
python -m pytest tests/ -v

# Test individual analyzers
python -m pytest tests/test_analyzer.py -v        # MHC-I tests
python -m pytest tests/test_mhc_ii_analyzer.py -v # MHC-II tests

# Benchmark performance
python tests/benchmark_analyzer.py
```

## Repository Structure

```
phlatcr_splicer/
├── phlatcr_splicer/              # Main package
│   ├── __init__.py               # Package initialization
│   ├── mhc_i_analyzer.py         # MHC-I analyzer
│   └── mhc_ii_analyzer.py        # MHC-II analyzer
├── scripts/                      # Command-line scripts
│   ├── __init__.py
│   └── main.py                   # Unified main entry point
├── tests/                        # Test suite
│   ├── test_analyzer.py          # MHC-I tests
│   ├── test_mhc_ii_analyzer.py   # MHC-II tests
│   └── benchmark_analyzer.py     # Performance tests
├── test_data/                    # Example PDB files
│   ├── 1oga.pdb                  # MHC-I example
│   ├── 4z7u.pdb                  # MHC-II example
│   └── 8ye4.pdb                  # Multi-complex example
├── examples/                     # Usage examples
├── docs/                         # Documentation
├── setup.py                      # Package installation
├── requirements.txt              # Dependencies
└── README.md                     # This file
```

## Algorithm Features

### Enhanced Complex Detection
- **Peptide-count based estimation** for reliable complex counting
- **Spatial clustering** for structures with many chains
- **Sophisticated pairing logic** understanding MHC-II α/β heterodimers

### Advanced Pattern Recognition  
- **Sequence motifs**: Highly specific patterns for each chain type
- **Length profiles**: Characteristic size ranges with optimal zones
- **Composition analysis**: Amino acid content scoring
- **Process-of-elimination**: Fallback classification for difficult sequences

### Multi-Complex Support
- **Automatic detection** of multiple complexes in single PDB files
- **Chain grouping** by spatial proximity and functional relationships
- **Complex numbering** for clear result organization

## Installation & Setup

### Dependencies

- Python 3.7+
- BioPython >= 1.79
- NumPy >= 1.19.0

### From Source

```bash
git clone https://github.com/yourusername/phlatcr_splicer.git
cd phlatcr_splicer
pip install -r requirements.txt
pip install -e .
```

### Verify Installation

```bash
# Test package installation
python -c "from phlatcr_splicer import pMHCITCRAnalyzer, pMHCIITCRAnalyzer; print('Package installed successfully!')"

# Test command-line tools
python scripts/main.py --help
phlatcr-analyze --help
mhc-ii-analyze --help
```

## Contributing

Please read [CONTRIBUTING.md](docs/CONTRIBUTING.md) for contribution guidelines.

## License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.

## Acknowledgments

- Xi'an Jiaotong-Liverpool University (XJTLU) Summer Undergraduate Research Fellowship (SURF)
- SURF Codes: 1335, 1415
