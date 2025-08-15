# pHLA-TCR Splicer: Unified TCR-pMHC Complex Analyzer

A sophisticated Python tool for analyzing TCR-pMHC protein complex structures from PDB and CIF files. This unified analyzer automatically handles both **MHC Class I** and **MHC Class II** complexes using advanced sequence alignment, spatial clustering, and pattern recognition techniques.

## 🚀 Key Features

- **Unified Analysis**: Single analyzer handles both MHC-I and MHC-II complexes automatically
- **Sequence Alignment**: High-confidence TCR α/β discrimination using reference sequences
- **Spatial Clustering**: DBSCAN-based detection of multiple complexes in single files
- **Multi-Format Support**: Works with both PDB and CIF file formats
- **Proximity-Based Reclassification**: Intelligent inference of unknown chains
- **Configurable Parameters**: 12+ tunable parameters for fine-grained control
- **Batch Processing**: Analyze multiple files with comprehensive CSV summaries
- **Chain Pairing Logic**: Sophisticated pairing of TCR dimers and pMHC complexes

## 📦 Installation

### Requirements
- Python 3.8+
- BioPython >= 1.79
- NumPy >= 1.21.0
- scikit-learn >= 1.0.0

### Install from Source

```bash
git clone https://github.com/skblnw/phlatcr_splicer.git
cd phlatcr_splicer
pip install -r requirements.txt
pip install -e .
```

### Verify Installation

```bash
# Test package import
python -c "from phlatcr_splicer import TCRpMHCAnalyzer; print('✅ Package installed successfully!')"

# Test CLI
python scripts/main.py --help
```

## 🎯 Quick Start

### Basic Usage

```bash
# Analyze a single PDB file
python scripts/main.py structure.pdb

# Analyze a CIF file with verbose output
python scripts/main.py structure.cif --verbose

# Batch process multiple files
python scripts/main.py *.pdb --batch-summary

# Save results to file
python scripts/main.py complex.pdb --output results.txt
```

### Advanced Usage with Parameters

```bash
# Tighter clustering for closely packed structures
python scripts/main.py structure.pdb --eps 45.0

# Stricter TCR alignment scoring
python scripts/main.py structure.pdb --align-score 25

# Custom distance thresholds
python scripts/main.py structure.pdb \
    --tcr-pair-dist 40.0 \
    --mhc1-pair-dist 45.0 \
    --pep-mhc1-dist 35.0
```

## 🧬 Chain Types Identified

The analyzer identifies the following chain types:

### Common to Both Complex Types
- **`TCR_ALPHA`**: T-cell receptor α chain (90-320 residues)
- **`TCR_BETA`**: T-cell receptor β chain (90-320 residues)
- **`PEPTIDE`**: Presented antigen (5-40 residues)

### MHC Class I Specific
- **`MHC_I_ALPHA`**: MHC-I heavy chain (240-300 residues)
- **`B2M`**: β2-microglobulin (80-120 residues)

### MHC Class II Specific
- **`MHC_II_ALPHA`**: MHC-II α chain (150-220 residues)
- **`MHC_II_BETA`**: MHC-II β chain (150-230 residues)

## 🔬 Algorithm Overview

### Three-Tier Chain Identification

1. **Length-Based Peptide Detection**: Peptides identified by characteristic length (5-40 residues)
2. **Alignment-Based TCR Classification**: TCR α/β chains discriminated using sequence alignment against reference constant regions
3. **Pattern-Based MHC Classification**: MHC and other chains identified through weighted pattern matching

### Spatial Organization

- **DBSCAN Clustering**: Groups chains into molecular complexes based on geometric centers
- **MHC-II Pairing Enforcement**: Ensures α/β heterodimers are kept together
- **Hierarchical Complex Assembly**: TCR dimers → pMHC complexes → full TCR-pMHC pairs

### Intelligent Reclassification

Proximity-based rules for unknown chains:
- UNKNOWN near MHC-I α → B2M
- UNKNOWN near B2M → MHC-I α
- UNKNOWN near MHC-II β → MHC-II α
- UNKNOWN near MHC-II α → MHC-II β

## 💻 Python API

```python
from phlatcr_splicer import TCRpMHCAnalyzer

# Initialize analyzer with default parameters
analyzer = TCRpMHCAnalyzer(verbose=True)

# Or with custom parameters
analyzer = TCRpMHCAnalyzer(
    eps=45.0,              # Tighter clustering
    align_score=25,        # Stricter alignment
    tcr_pair_dist=40.0,    # Closer TCR pairing
    verbose=True
)

# Analyze a structure
results, chain_map = analyzer.analyze_pdb("structure.pdb")

# Process results
for complex_name, complex_data in results.items():
    chains = complex_data['chains']
    pairs = complex_data['pairs']
    
    print(f"{complex_name}:")
    for chain_id, chain_type in chains.items():
        print(f"  Chain {chain_id}: {chain_type}")
    
    for pair in pairs:
        print(f"  Complex type: {pair['type']}")
```

## 📊 Output Format

### Console Output Example

```
📁 Found 1 file(s) to analyze
🧬 Analyzing structure.pdb...
============================================================

📊 Analysis Results for structure.pdb:
--------------------------------------------------

Complex 1:
  Chain Assignments:
    Chain A: MHC-I Heavy Chain
    Chain B: β2-microglobulin
    Chain C: Peptide
    Chain D: TCR-α
    Chain E: TCR-β
  Identified Complexes:
    Complex #1: TCR-α[D] + TCR-β[E] :: MHC-I[A] + B2M[B] + Pep[C]

✅ Summary:
  Total chains identified: 5
  Total complexes found: 1
  MHC-I complexes: 1

✨ Analysis complete! Processed 1 file(s)
```

### CSV Output Format

For batch processing with `--batch-summary`:

| PDB File | Complex ID | Complex Type | TCR_Alpha | TCR_Beta | MHC_Alpha | MHC_Beta | B2M | Peptide |
|----------|------------|--------------|-----------|----------|-----------|----------|-----|---------|
| 1oga.pdb | complex_1_1 | MHC-I | D | E | A |  | B | C |
| 4z7u.pdb | complex_1_1 | MHC-II | E | F | A | B |  | I |

## ⚙️ Configuration Parameters

| Parameter | Default | Description |
|-----------|---------|-------------|
| `--eps` | 60.0 | DBSCAN clustering distance (Å) |
| `--align-score` | 20 | Minimum alignment score for TCR classification |
| `--align-ratio` | 1.5 | Score ratio for TCR α/β discrimination |
| `--tcr-pair-dist` | 50.0 | Max distance for TCR α/β pairing (Å) |
| `--mhc1-pair-dist` | 50.0 | Max distance for MHC-I/B2M pairing (Å) |
| `--mhc2-pair-dist` | 50.0 | Max distance for MHC-II α/β pairing (Å) |
| `--pep-mhc1-dist` | 40.0 | Max distance for Peptide/MHC-I pairing (Å) |
| `--pep-mhc2-dist` | 60.0 | Max distance for Peptide/MHC-II pairing (Å) |
| `--tcr-pmhc-dist` | 150.0 | Max distance for TCR/pMHC pairing (Å) |
| `--tcr-pmhc-dist-unpaired` | 120.0 | Max distance for single-chain TCR/pMHC (Å) |
| `--reclassify-dist` | 50.0 | Max distance for reclassifying UNKNOWN chains (Å) |

## 🧪 Testing

```bash
# Run all tests
python -m pytest tests/ -v

# Run specific test modules
python -m pytest tests/test_unified_analyzer.py -v
python -m pytest tests/test_cli.py -v
python -m pytest tests/test_integration.py -v

# Run with coverage
python -m pytest tests/ --cov=phlatcr_splicer --cov-report=html
```

## 📁 Repository Structure

```
phlatcr_splicer/
├── phlatcr_splicer/
│   ├── __init__.py           # Package initialization
│   └── unified_analyzer.py   # Main analyzer module
├── scripts/
│   └── main.py              # CLI entry point
├── tests/
│   ├── test_unified_analyzer.py  # Unit tests
│   ├── test_cli.py              # CLI tests
│   └── test_integration.py      # Integration tests
├── examples/
│   ├── basic_usage.py           # Simple examples
│   ├── batch_processing.py      # Batch analysis
│   └── advanced_parameters.py   # Parameter tuning
├── docs/
│   ├── API.md                   # API documentation
│   ├── ALGORITHM.md             # Algorithm details
│   └── MIGRATION.md             # Migration guide
├── test_data/                   # Example PDB/CIF files
├── setup.py                     # Package configuration
├── requirements.txt             # Dependencies
├── CHANGELOG.md                 # Version history
└── README.md                    # This file
```

## 🔄 Migration from Previous Versions

If you're upgrading from the dual-analyzer system (v0.x), see [MIGRATION.md](docs/MIGRATION.md) for details.

Key changes:
- No more `--type` or `--auto` flags needed
- Single `TCRpMHCAnalyzer` class replaces `pMHCITCRAnalyzer` and `pMHCIITCRAnalyzer`
- New parameters for fine-tuning analysis
- CIF file support added

## 📝 Citation

If you use this tool in your research, please cite:

```bibtex
@software{phlatcr_splicer,
  author = {skblnw},
  title = {pHLA-TCR Splicer: Unified TCR-pMHC Complex Analyzer},
  year = {2024},
  url = {https://github.com/skblnw/phlatcr_splicer}
}
```

## 🤝 Contributing

Contributions are welcome! Please feel free to submit issues or pull requests.

## 📄 License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.

## 🙏 Acknowledgments

- Xi'an Jiaotong-Liverpool University (XJTLU) Summer Undergraduate Research Fellowship (SURF)
- SURF Codes: 1335, 1415

## 📧 Contact

- Author: skblnw
- Email: skblnw@github.com
- GitHub: [@skblnw](https://github.com/skblnw)