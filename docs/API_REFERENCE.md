# API Reference

## pMHCITCRAnalyzer Class

The main class for analyzing pHLA-TCR complex structures.

### Constructor

```python
pMHCITCRAnalyzer(verbose: bool = True)
```

**Parameters:**
- `verbose` (bool): Whether to print verbose output during analysis

### Methods

#### analyze_pdb(pdb_file: str) -> Dict[str, str]

Analyze a PDB file and identify chain types.

**Parameters:**
- `pdb_file` (str): Path to the PDB file

**Returns:**
- `Dict[str, str]`: Dictionary mapping chain IDs to chain types

**Raises:**
- `FileNotFoundError`: If the PDB file doesn't exist
- `Exception`: For various PDB parsing or analysis errors

**Example:**
```python
analyzer = pMHCITCRAnalyzer()
result = analyzer.analyze_pdb("complex.pdb")
# Returns: {'A': 'tcr_alpha', 'B': 'tcr_beta', 'H': 'mhc_heavy', 'L': 'b2m', 'P': 'peptide'}
```

#### save_analysis_report(assignments: Dict[str, str], output_file: str, chain_info: Dict = None)

Save detailed analysis report to file.

**Parameters:**
- `assignments` (Dict[str, str]): Chain assignments from analyze_pdb()
- `output_file` (str): Output file path
- `chain_info` (Dict, optional): Additional chain information

**Example:**
```python
analyzer.save_analysis_report(result, "report.txt")
```

### Chain Types

The analyzer identifies the following chain types:

- **`peptide`**: Peptide antigen (typically 8-15 residues)
- **`mhc_heavy`**: MHC class I heavy chain (270-380 residues)
- **`b2m`**: β2-microglobulin (95-105 residues)  
- **`tcr_alpha`**: TCR alpha chain (200-250 residues)
- **`tcr_beta`**: TCR beta chain (240-290 residues)
- **`unknown`**: Unidentified or ambiguous chains

### Internal Methods

#### Scoring Functions

These methods calculate confidence scores for each chain type:

- `_score_peptide(sequence, length, properties) -> float`
- `_score_b2m(sequence, length, properties) -> float`
- `_score_mhc(sequence, length, properties) -> float`
- `_score_tcr_alpha(sequence, length) -> float`
- `_score_tcr_beta(sequence, length) -> float`

#### Pattern Loading

- `_load_mhc_patterns() -> Dict`: Load MHC heavy chain patterns
- `_load_tcr_patterns() -> Dict`: Load TCR alpha/beta patterns  
- `_load_b2m_patterns() -> Dict`: Load β2-microglobulin patterns

#### Utility Functions

- `_extract_chain_info(structure) -> Dict`: Extract chain information from PDB structure
- `_get_chain_sequence(chain) -> str`: Get amino acid sequence from chain
- `_three_to_one(three_letter) -> str`: Convert 3-letter to 1-letter amino acid codes
- `_calculate_chain_properties(sequence, chain) -> Dict`: Calculate molecular properties
- `_classify_chains(chain_info) -> Dict[str, str]`: Main classification logic
- `_determine_chain_type(sequence, length, properties, used_types) -> str`: Determine single chain type
- `_validate_assignments(assignments, chain_info) -> Dict[str, str]`: Validate and correct assignments

## Command Line Interface

```bash
python phlatcr_analyzer.py [options] pdb_file
```

**Arguments:**
- `pdb_file`: Input PDB file path

**Options:**
- `-o, --output`: Output report file
- `-v, --verbose`: Verbose output

**Examples:**
```bash
# Basic analysis
python phlatcr_analyzer.py complex.pdb

# With verbose output and report
python phlatcr_analyzer.py -v -o report.txt complex.pdb
```

## Pattern Definitions

### MHC Heavy Chain Patterns

```python
{
    'length_range': (270, 380),
    'conserved_motifs': ['GSHSMRY', 'WSDRVI', 'FKAFLKQ', 'TGAASC', 'RGEC'],
    'domain_boundaries': {
        'alpha1': (1, 90),
        'alpha2': (91, 180), 
        'alpha3': (181, 276)
    },
    'disulfide_cysteines': [25, 101, 164, 236]
}
```

### TCR Patterns

```python
{
    'alpha': {
        'length_range': (200, 250),
        'conserved_motifs': ['FGXGT', 'WYQQKP', 'YKFK', 'TLTIS'],
        'cdr3_patterns': ['CAS', 'CAV', 'CAI']
    },
    'beta': {
        'length_range': (240, 290),
        'conserved_motifs': ['FGGGT', 'WYQQKP', 'MGIGV', 'SVGD'],
        'cdr3_patterns': ['CAS', 'CAW', 'CAT']
    }
}
```

### β2-Microglobulin Patterns

```python
{
    'length_range': (95, 105),
    'conserved_motifs': ['IQRTPK', 'VNHVTL', 'SRNLTKDR'],
    'molecular_weight': (11000, 13000)
}
```

## Error Handling

The analyzer handles several error conditions:

- **File not found**: Raises `FileNotFoundError` for missing PDB files
- **Invalid PDB format**: BioPython parsing exceptions are caught and reported
- **Empty chains**: Chains with no standard amino acids are handled gracefully
- **Ambiguous assignments**: Validation logic prevents duplicate unique chain types

## Performance Considerations

- **Memory usage**: ~100MB for typical complexes
- **Speed**: <1 second analysis time for most complexes
- **Scalability**: Linear scaling with number of chains

## Extending the Analyzer

### Adding New Chain Types

1. Add patterns to pattern loading methods
2. Create new scoring function
3. Update classification logic
4. Add validation rules
5. Include in tests

### Custom Pattern Matching

Subclass `pMHCITCRAnalyzer` and override pattern loading methods:

```python
class CustomAnalyzer(pMHCITCRAnalyzer):
    def _load_custom_patterns(self):
        return {
            'length_range': (100, 200),
            'conserved_motifs': ['CUSTOM', 'MOTIF'],
        }
```