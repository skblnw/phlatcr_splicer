# pHLA-TCR Complex Structure Analyzer - User Guide

## Table of Contents

1. [Introduction](#introduction)
2. [Installation](#installation)
3. [Quick Start](#quick-start)
4. [Command Line Usage](#command-line-usage)
5. [Python API Usage](#python-api-usage)
6. [Understanding Results](#understanding-results)
7. [Troubleshooting](#troubleshooting)
8. [Advanced Usage](#advanced-usage)
9. [Examples](#examples)

## Introduction

The pHLA-TCR Complex Structure Analyzer is a Python tool designed to automatically identify and classify protein chains in pHLA-TCR (peptide-HLA-T Cell Receptor) complexes from PDB structure files. It's particularly useful when dealing with structures where chain IDs are inconsistent or non-standard.

### What it identifies:
- **Peptide antigen**: The short peptide presented by MHC
- **MHC heavy chain**: The major histocompatibility complex class I heavy chain
- **β2-microglobulin (b2m)**: The MHC light chain
- **TCR alpha chain**: T cell receptor alpha chain
- **TCR beta chain**: T cell receptor beta chain

## Installation

### Requirements
- Python 3.7 or higher
- BioPython ≥ 1.79
- NumPy ≥ 1.21.0
- SciPy ≥ 1.7.0

### Install from source
```bash
git clone https://github.com/yourusername/phlatcr_splicer.git
cd phlatcr_splicer
pip install -r requirements.txt
```

### Verify installation
```bash
python phlatcr_analyzer.py --help
```

## Quick Start

### Basic command line usage:
```bash
python phlatcr_analyzer.py your_complex.pdb
```

### Basic Python usage:
```python
from phlatcr_analyzer import pHLATCRAnalyzer

analyzer = pHLATCRAnalyzer()
result = analyzer.analyze_pdb("your_complex.pdb")
print(result)
```

## Command Line Usage

### Basic Analysis
```bash
python phlatcr_analyzer.py complex.pdb
```

### With Verbose Output
```bash
python phlatcr_analyzer.py -v complex.pdb
```

### Save Report
```bash
python phlatcr_analyzer.py -o analysis_report.txt complex.pdb
```

### Full Command Options
```bash
python phlatcr_analyzer.py [-h] [-o OUTPUT] [-v] pdb_file

Arguments:
  pdb_file              Input PDB file

Options:
  -h, --help           Show help message
  -o, --output OUTPUT  Output report file
  -v, --verbose        Verbose output
```

## Python API Usage

### Basic Analysis
```python
from phlatcr_analyzer import pHLATCRAnalyzer

# Initialize analyzer
analyzer = pHLATCRAnalyzer(verbose=True)

# Analyze PDB file
assignments = analyzer.analyze_pdb("complex.pdb")

# Print results
for chain_id, chain_type in assignments.items():
    print(f"Chain {chain_id}: {chain_type}")
```

### Batch Processing
```python
import os
from phlatcr_analyzer import pHLATCRAnalyzer

analyzer = pHLATCRAnalyzer(verbose=False)
results = {}

# Process multiple files
pdb_files = ["complex1.pdb", "complex2.pdb", "complex3.pdb"]

for pdb_file in pdb_files:
    if os.path.exists(pdb_file):
        try:
            result = analyzer.analyze_pdb(pdb_file)
            results[pdb_file] = result
            print(f"✓ Processed {pdb_file}")
        except Exception as e:
            print(f"✗ Error processing {pdb_file}: {e}")

# Save batch report
with open("batch_results.txt", "w") as f:
    for pdb_file, result in results.items():
        f.write(f"{pdb_file}:\n")
        for chain_id, chain_type in result.items():
            f.write(f"  {chain_id}: {chain_type}\n")
        f.write("\n")
```

### Custom Analysis
```python
class DetailedAnalyzer(pHLATCRAnalyzer):
    def analyze_with_confidence(self, pdb_file):
        # Get basic assignments
        assignments = self.analyze_pdb(pdb_file)
        
        # Add confidence scores
        structure = self.parser.get_structure('complex', pdb_file)
        chain_info = self._extract_chain_info(structure)
        
        detailed_results = {}
        for chain_id, chain_type in assignments.items():
            info = chain_info[chain_id]
            
            # Calculate confidence score
            if chain_type == 'peptide':
                confidence = self._score_peptide(info['sequence'], info['length'], info['properties'])
            elif chain_type == 'b2m':
                confidence = self._score_b2m(info['sequence'], info['length'], info['properties'])
            # ... etc for other types
            else:
                confidence = 0.0
            
            detailed_results[chain_id] = {
                'type': chain_type,
                'confidence': confidence,
                'length': info['length']
            }
        
        return detailed_results

# Use custom analyzer
analyzer = DetailedAnalyzer()
detailed_result = analyzer.analyze_with_confidence("complex.pdb")
```

## Understanding Results

### Output Format
The analyzer returns a dictionary mapping chain IDs to chain types:

```python
{
    'A': 'tcr_alpha',
    'B': 'tcr_beta', 
    'H': 'mhc_heavy',
    'L': 'b2m',
    'P': 'peptide'
}
```

### Chain Types
- **`peptide`**: Short peptide antigen (typically 8-15 residues)
- **`mhc_heavy`**: MHC class I heavy chain (270-380 residues)
- **`b2m`**: β2-microglobulin (95-105 residues)
- **`tcr_alpha`**: TCR alpha chain (200-250 residues)
- **`tcr_beta`**: TCR beta chain (240-290 residues)
- **`unknown`**: Could not be confidently assigned

### Confidence Levels
The analyzer uses multiple criteria for chain identification:
- **High confidence**: Multiple pattern matches, correct length, good score
- **Medium confidence**: Some pattern matches, reasonable length
- **Low confidence**: Length-based assignment only, falls back to 'unknown'

## Troubleshooting

### Common Issues

#### 1. "FileNotFoundError: PDB file not found"
- **Cause**: Incorrect file path
- **Solution**: Check file path and ensure file exists

#### 2. "Invalid or missing coordinate(s)"
- **Cause**: Corrupted or incorrectly formatted PDB file
- **Solution**: Verify PDB file format, try re-downloading structure

#### 3. "No chains found" or empty results
- **Cause**: PDB file contains no standard amino acid residues
- **Solution**: Check if file contains protein chains vs. other molecules

#### 4. Multiple chains assigned same type
- **Cause**: Ambiguous structure or similar chain lengths
- **Solution**: The analyzer includes validation to resolve this automatically

#### 5. Poor chain assignments
- **Cause**: Heavily modified or non-standard sequences
- **Solution**: 
  - Use verbose mode to see scoring details
  - Consider manual review for engineered constructs
  - Check if structure contains complete chains

### Debug Mode
```python
# Enable verbose output for debugging
analyzer = pHLATCRAnalyzer(verbose=True)
result = analyzer.analyze_pdb("problematic.pdb")
```

### Validating Results
```python
# Check assignments make biological sense
def validate_phla_tcr_complex(assignments):
    expected_types = {'peptide', 'mhc_heavy', 'b2m', 'tcr_alpha', 'tcr_beta'}
    found_types = set(assignments.values())
    
    missing = expected_types - found_types
    unexpected = found_types - expected_types - {'unknown'}
    
    print(f"Missing chain types: {missing}")
    print(f"Unexpected chain types: {unexpected}")
    
    return len(missing) == 0 and len(unexpected) == 0

# Usage
is_valid = validate_phla_tcr_complex(result)
```

## Advanced Usage

### Custom Pattern Matching
```python
class CustomAnalyzer(pHLATCRAnalyzer):
    def _load_custom_patterns(self):
        """Add patterns for custom chain types."""
        return {
            'my_protein': {
                'length_range': (150, 200),
                'conserved_motifs': ['SPECIAL', 'MOTIF'],
                'molecular_weight': (15000, 20000)
            }
        }
    
    def _score_my_protein(self, sequence, length, properties):
        """Custom scoring function."""
        score = 0.0
        if 150 <= length <= 200:
            score += 0.5
        if 'SPECIAL' in sequence:
            score += 0.3
        return min(score, 1.0)
```

### Integration with Structure Analysis
```python
from Bio.PDB import *

def analyze_contacts(pdb_file, assignments):
    """Find inter-chain contacts."""
    parser = PDBParser()
    structure = parser.get_structure('complex', pdb_file)
    
    # Get peptide and MHC chains
    peptide_chain = None
    mhc_chain = None
    
    for chain in structure.get_chains():
        chain_id = chain.get_id()
        if assignments.get(chain_id) == 'peptide':
            peptide_chain = chain
        elif assignments.get(chain_id) == 'mhc_heavy':
            mhc_chain = chain
    
    if peptide_chain and mhc_chain:
        # Analyze peptide-MHC contacts
        contacts = find_contacts(peptide_chain, mhc_chain, cutoff=4.0)
        return contacts
    
    return []
```

## Examples

### Example 1: Single Complex Analysis
```python
from phlatcr_analyzer import pHLATCRAnalyzer

# Analyze a single pHLA-TCR complex
analyzer = pHLATCRAnalyzer(verbose=True)
result = analyzer.analyze_pdb("1ao7.pdb")

print("Chain assignments:")
for chain_id, chain_type in result.items():
    print(f"  Chain {chain_id}: {chain_type}")

# Save detailed report
analyzer.save_analysis_report(result, "1ao7_analysis.txt")
```

### Example 2: Batch Processing with Error Handling
```python
import glob
from phlatcr_analyzer import pHLATCRAnalyzer

analyzer = pHLATCRAnalyzer(verbose=False)

# Process all PDB files in directory
pdb_files = glob.glob("*.pdb")
successful = 0
failed = 0

for pdb_file in pdb_files:
    try:
        result = analyzer.analyze_pdb(pdb_file)
        
        # Check if we found a complete complex
        chain_types = set(result.values())
        expected = {'peptide', 'mhc_heavy', 'b2m', 'tcr_alpha', 'tcr_beta'}
        
        if expected.issubset(chain_types):
            print(f"✓ Complete complex: {pdb_file}")
            successful += 1
        else:
            print(f"⚠ Incomplete complex: {pdb_file}")
            print(f"   Found: {chain_types}")
            
    except Exception as e:
        print(f"✗ Failed: {pdb_file} - {e}")
        failed += 1

print(f"\nSummary: {successful} successful, {failed} failed")
```

### Example 3: Integration with Data Analysis
```python
import pandas as pd
from phlatcr_analyzer import pHLATCRAnalyzer

def analyze_complex_dataset(pdb_files):
    """Analyze a dataset of pHLA-TCR complexes."""
    analyzer = pHLATCRAnalyzer(verbose=False)
    
    results = []
    for pdb_file in pdb_files:
        try:
            assignments = analyzer.analyze_pdb(pdb_file)
            
            # Extract chain lengths and types
            structure = analyzer.parser.get_structure('complex', pdb_file)
            chain_info = analyzer._extract_chain_info(structure)
            
            row = {'pdb_file': pdb_file}
            for chain_id, chain_type in assignments.items():
                row[f'{chain_type}_length'] = chain_info[chain_id]['length']
                row[f'{chain_type}_chain'] = chain_id
            
            results.append(row)
            
        except Exception as e:
            print(f"Error with {pdb_file}: {e}")
    
    # Convert to DataFrame for analysis
    df = pd.DataFrame(results)
    return df

# Usage
pdb_files = ["complex1.pdb", "complex2.pdb", "complex3.pdb"]
df = analyze_complex_dataset(pdb_files)
print(df.describe())
```

For more examples, see `example_usage.py` in the repository.