# Changelog

All notable changes to the pHLA-TCR Splicer project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [1.0.0] - 2024-01-XX

### 🚀 Major Release - Unified Analyzer

This release introduces a completely redesigned unified analyzer that handles both MHC-I and MHC-II complexes in a single system, with significant improvements in accuracy and functionality.

### Added
- **Unified TCRpMHCAnalyzer**: Single analyzer class for both MHC-I and MHC-II complexes
- **Sequence Alignment**: High-confidence TCR α/β discrimination using reference sequences
- **DBSCAN Clustering**: Automatic detection of multiple complexes in single files
- **CIF Support**: Full support for mmCIF format files with chain ID mapping
- **Proximity-Based Reclassification**: Intelligent inference of unknown chains based on spatial relationships
- **12+ Configuration Parameters**: Fine-grained control over analysis behavior
- **Batch Processing**: Analyze multiple files with CSV summary output
- **Comprehensive Test Suite**: New tests for CLI, integration, and performance
- **Enhanced Examples**: Basic usage, batch processing, and parameter tuning examples

### Changed
- **Single Entry Point**: Removed need for `--type` and `--auto` flags
- **Improved Accuracy**: Three-tier chain identification strategy
- **Better Performance**: Optimized clustering and pairing algorithms
- **Cleaner API**: Simplified to single `TCRpMHCAnalyzer` class
- **Updated Dependencies**: Added scikit-learn for clustering

### Removed
- **Deprecated Analyzers**: Removed separate `pMHCITCRAnalyzer` and `pMHCIITCRAnalyzer` classes
- **Auto-Detection Logic**: No longer needed with unified approach
- **Redundant Tests**: Removed old test files for deprecated analyzers
- **Dual Analyzer Examples**: Removed obsolete dual analyzer example

### Fixed
- **Multi-Complex Detection**: Improved handling of multiple complexes in single files
- **MHC-II Pairing**: Better enforcement of α/β heterodimer relationships
- **Unknown Chain Handling**: More robust reclassification system

### Migration
See [MIGRATION.md](docs/MIGRATION.md) for detailed migration instructions from v0.x to v1.0.0.

## [0.2.0] - 2023-XX-XX

### Added
- Initial unified analyzer prototype
- Deprecation warnings for old analyzers

## [0.1.0] - 2023-XX-XX

### Initial Release
- Separate MHC-I and MHC-II analyzers
- Basic pattern recognition
- Command-line interface
- Auto-detection between MHC types

---

For more details on each release, see the [GitHub releases page](https://github.com/skblnw/phlatcr_splicer/releases).