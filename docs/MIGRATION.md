# Migration Guide: v0.x to v1.0.0

This guide helps you migrate from the dual-analyzer system (v0.x) to the unified analyzer (v1.0.0).

## Overview of Changes

The v1.0.0 release introduces a unified analyzer that handles both MHC-I and MHC-II complexes automatically. This eliminates the need for separate analyzers and manual type detection.

## Quick Migration

### Command Line

**Old way (v0.x):**
```bash
# Had to specify type
python scripts/main.py --type mhc-i input.pdb
python scripts/main.py --type mhc-ii input.pdb

# Or use auto-detection
python scripts/main.py --auto input.pdb
```

**New way (v1.0.0):**
```bash
# Automatic detection - no type needed!
python scripts/main.py input.pdb
```

### Python API

**Old way (v0.x):**
```python
from phlatcr_splicer import pMHCITCRAnalyzer, pMHCIITCRAnalyzer

# Had to choose analyzer
mhc_i_analyzer = pMHCITCRAnalyzer(verbose=True)
mhc_ii_analyzer = pMHCIITCRAnalyzer(verbose=True)

# Analyze based on type
if complex_type == "mhc-i":
    results = mhc_i_analyzer.analyze_pdb("file.pdb")
else:
    results = mhc_ii_analyzer.analyze_pdb("file.pdb")
```

**New way (v1.0.0):**
```python
from phlatcr_splicer import TCRpMHCAnalyzer

# Single analyzer for everything
analyzer = TCRpMHCAnalyzer(verbose=True)
results, chain_map = analyzer.analyze_pdb("file.pdb")
```

## Detailed Changes

### 1. Class Names

| Old (v0.x) | New (v1.0.0) |
|------------|--------------|
| `pMHCITCRAnalyzer` | `TCRpMHCAnalyzer` |
| `pMHCIITCRAnalyzer` | `TCRpMHCAnalyzer` |

### 2. Return Values

**Old way:** Single dictionary of chain assignments
```python
results = analyzer.analyze_pdb("file.pdb")
# Returns: {'A': 'mhc_heavy', 'B': 'b2m', ...}
```

**New way:** Tuple of (results, chain_map)
```python
results, chain_map = analyzer.analyze_pdb("file.pdb")
# results: Complex groupings with chains and pairs
# chain_map: Mapping for CIF files (empty dict for PDB)
```

### 3. Result Structure

**Old structure:**
```python
{
    'A': 'mhc_heavy',
    'B': 'b2m',
    'C': 'peptide',
    'D': 'tcr_alpha',
    'E': 'tcr_beta'
}
```

**New structure:**
```python
{
    'complex_1': {
        'chains': {
            'A': 'MHC_I_ALPHA',
            'B': 'B2M',
            'C': 'PEPTIDE',
            'D': 'TCR_ALPHA',
            'E': 'TCR_BETA'
        },
        'pairs': [
            {
                'type': 'MHC-I',
                'mhc_alpha': 'A',
                'b2m': 'B',
                'peptide': 'C',
                'tcr_alpha': 'D',
                'tcr_beta': 'E'
            }
        ]
    }
}
```

### 4. Chain Type Names

| Old (v0.x) | New (v1.0.0) |
|------------|--------------|
| `mhc_heavy` | `MHC_I_ALPHA` |
| `b2m` | `B2M` |
| `peptide` | `PEPTIDE` |
| `tcr_alpha` | `TCR_ALPHA` |
| `tcr_beta` | `TCR_BETA` |
| `mhc_ii_alpha` | `MHC_II_ALPHA` |
| `mhc_ii_beta` | `MHC_II_BETA` |

### 5. New Parameters

The unified analyzer adds many tunable parameters:

```python
analyzer = TCRpMHCAnalyzer(
    eps=45.0,              # DBSCAN clustering distance
    align_score=25,        # TCR alignment threshold
    tcr_pair_dist=40.0,    # TCR α/β pairing distance
    # ... and more
)
```

## Step-by-Step Migration

### Step 1: Update Imports

```python
# Remove old imports
# from phlatcr_splicer import pMHCITCRAnalyzer, pMHCIITCRAnalyzer

# Add new import
from phlatcr_splicer import TCRpMHCAnalyzer
```

### Step 2: Update Analyzer Initialization

```python
# Replace old analyzers
# analyzer = pMHCITCRAnalyzer(verbose=True)

# With unified analyzer
analyzer = TCRpMHCAnalyzer(verbose=True)
```

### Step 3: Update Result Processing

```python
# Old way
# results = analyzer.analyze_pdb("file.pdb")
# for chain_id, chain_type in results.items():
#     print(f"Chain {chain_id}: {chain_type}")

# New way
results, chain_map = analyzer.analyze_pdb("file.pdb")
for complex_name, complex_data in results.items():
    chains = complex_data['chains']
    pairs = complex_data['pairs']
    
    for chain_id, chain_type in chains.items():
        print(f"Chain {chain_id}: {chain_type}")
    
    for pair in pairs:
        print(f"Complex type: {pair['type']}")
```

### Step 4: Update Chain Type Checks

```python
# Old way
# if chain_type == 'mhc_heavy':
#     process_mhc_heavy_chain()

# New way
if chain_type == 'MHC_I_ALPHA':
    process_mhc_heavy_chain()
```

### Step 5: Remove Type Detection

Remove any code that detects or switches between MHC-I and MHC-II:

```python
# Remove this kind of logic
# if has_b2m:
#     analyzer = pMHCITCRAnalyzer()
# else:
#     analyzer = pMHCIITCRAnalyzer()

# Just use unified analyzer
analyzer = TCRpMHCAnalyzer()
```

## New Features to Explore

### 1. Multi-Complex Support

The new analyzer automatically detects multiple complexes:

```python
results, _ = analyzer.analyze_pdb("multi_complex.pdb")
print(f"Found {len(results)} complexes")
```

### 2. CIF File Support

Now supports CIF files with chain mapping:

```python
results, chain_map = analyzer.analyze_pdb("structure.cif")
if chain_map:
    print("Chain ID mapping:", chain_map)
```

### 3. Fine-Tuning Parameters

Experiment with parameters for better results:

```python
# For tightly packed structures
analyzer = TCRpMHCAnalyzer(eps=45.0)

# For stricter TCR identification
analyzer = TCRpMHCAnalyzer(align_score=25)
```

## Troubleshooting

### Issue: ImportError for old classes

**Error:**
```python
ImportError: cannot import name 'pMHCITCRAnalyzer'
```

**Solution:** Update to use `TCRpMHCAnalyzer` instead.

### Issue: Results structure different

**Problem:** Code expects flat dictionary, gets nested structure.

**Solution:** Update result processing to handle new structure with complexes and pairs.

### Issue: Chain types not matching

**Problem:** String comparisons failing.

**Solution:** Update chain type names (e.g., `mhc_heavy` → `MHC_I_ALPHA`).

## Need Help?

- Check the [examples](../examples/) directory for updated usage patterns
- See the [API documentation](API.md) for detailed method signatures
- Open an issue on [GitHub](https://github.com/skblnw/phlatcr_splicer/issues)

## Benefits of Migration

1. **Simpler Code**: No need to manage two analyzers
2. **Better Accuracy**: Sequence alignment and clustering
3. **More Features**: CIF support, multi-complex detection
4. **Future-Proof**: Active development on unified analyzer
5. **Better Performance**: Optimized algorithms