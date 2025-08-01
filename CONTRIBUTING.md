# Contributing to pHLA-TCR Splicer

Thank you for your interest in contributing to the pHLA-TCR Complex Structure Analyzer!

## Development Setup

1. **Clone the repository**:
   ```bash
   git clone https://github.com/yourusername/phlatcr_splicer.git
   cd phlatcr_splicer
   ```

2. **Install dependencies**:
   ```bash
   pip install -r requirements.txt
   ```

3. **Run tests**:
   ```bash
   python test_analyzer.py
   python benchmark_analyzer.py
   ```

## Code Organization

- `phlatcr_analyzer.py` - Main analyzer class and algorithms
- `test_analyzer.py` - Unit tests and validation
- `benchmark_analyzer.py` - Performance benchmarking
- `example_usage.py` - Usage examples and demos

## Adding New Features

### Chain Type Detection

To add a new chain type:

1. Add patterns to the appropriate `_load_*_patterns()` method
2. Create a new scoring function `_score_new_type()`
3. Update `_determine_chain_type()` to include the new type
4. Add tests for the new functionality

### Pattern Matching

Pattern matching is based on:
- Sequence length ranges
- Conserved motifs
- Molecular weight
- Structural features (e.g., disulfide bonds)

### Testing

Always add tests for new features:
- Unit tests for individual functions
- Integration tests for complete workflows
- Benchmark tests for performance validation

## Code Style

- Follow PEP 8 style guidelines
- Use type hints where appropriate
- Write descriptive docstrings
- Keep functions focused and modular

## Pull Request Process

1. Fork the repository
2. Create a feature branch (`git checkout -b feature/new-feature`)
3. Make your changes
4. Add tests for new functionality
5. Ensure all tests pass
6. Submit a pull request with a clear description

## Reporting Issues

When reporting issues, please include:
- Python version
- BioPython version
- Sample PDB file (if applicable)
- Full error traceback
- Steps to reproduce

## Feature Requests

We welcome feature requests! Please describe:
- The biological use case
- Expected input/output
- Any relevant references

## Code of Conduct

Be respectful and constructive in all interactions. This project aims to be welcoming to contributors from all backgrounds.