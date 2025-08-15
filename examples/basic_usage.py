#!/usr/bin/env python3
"""
Basic usage examples for the Unified TCR-pMHC Analyzer

This script demonstrates simple usage patterns for analyzing
TCR-pMHC complex structures.

Author: skblnw
"""

from phlatcr_splicer import TCRpMHCAnalyzer


def analyze_single_file():
    """Basic analysis of a single PDB file."""
    print("=" * 60)
    print("Basic Single File Analysis")
    print("=" * 60)
    
    # Initialize analyzer with default parameters
    analyzer = TCRpMHCAnalyzer(verbose=True)
    
    # Analyze a PDB file
    pdb_file = "test_data/1oga.pdb"  # Replace with your file
    
    try:
        results, chain_map = analyzer.analyze_pdb(pdb_file)
        
        print(f"\nAnalysis complete for {pdb_file}")
        print(f"Found {len(results)} complex(es)\n")
        
        # Process results
        for complex_name, complex_data in results.items():
            print(f"{complex_name}:")
            
            # Show chain assignments
            chains = complex_data['chains']
            print("  Chain assignments:")
            for chain_id, chain_type in sorted(chains.items()):
                print(f"    Chain {chain_id}: {chain_type}")
            
            # Show paired complexes
            pairs = complex_data['pairs']
            print("  Paired complexes:")
            for i, pair in enumerate(pairs, 1):
                print(f"    Complex #{i}:")
                print(f"      Type: {pair.get('type', 'Unknown')}")
                if pair.get('tcr_alpha'):
                    print(f"      TCR-α: {pair['tcr_alpha']}")
                if pair.get('tcr_beta'):
                    print(f"      TCR-β: {pair['tcr_beta']}")
                if pair.get('mhc_alpha'):
                    print(f"      MHC-α: {pair['mhc_alpha']}")
                if pair.get('mhc_beta'):
                    print(f"      MHC-β: {pair['mhc_beta']}")
                if pair.get('b2m'):
                    print(f"      B2M: {pair['b2m']}")
                if pair.get('peptide'):
                    print(f"      Peptide: {pair['peptide']}")
            print()
    
    except FileNotFoundError:
        print(f"File {pdb_file} not found. Please provide a valid PDB file.")
    except Exception as e:
        print(f"Error analyzing file: {e}")


def analyze_with_custom_parameters():
    """Analysis with custom parameters for fine-tuning."""
    print("=" * 60)
    print("Analysis with Custom Parameters")
    print("=" * 60)
    
    # Initialize with custom parameters
    analyzer = TCRpMHCAnalyzer(
        eps=45.0,               # Tighter clustering for closely packed structures
        align_score=25,         # Stricter TCR alignment threshold
        tcr_pair_dist=40.0,     # Require TCR chains to be closer
        verbose=False           # Less output
    )
    
    pdb_file = "test_data/complex.pdb"  # Replace with your file
    
    try:
        results, chain_map = analyzer.analyze_pdb(pdb_file)
        
        print(f"\nAnalysis with custom parameters complete")
        print(f"Parameters used:")
        print(f"  - EPS (clustering): 45.0 Å")
        print(f"  - Alignment score: 25")
        print(f"  - TCR pair distance: 40.0 Å")
        print(f"\nResults: {len(results)} complex(es) found")
        
        # Summary statistics
        total_chains = sum(len(c['chains']) for c in results.values())
        total_pairs = sum(len(c['pairs']) for c in results.values())
        
        print(f"Total chains: {total_chains}")
        print(f"Total paired complexes: {total_pairs}")
    
    except FileNotFoundError:
        print(f"File {pdb_file} not found. Please provide a valid PDB file.")
    except Exception as e:
        print(f"Error: {e}")


def analyze_cif_file():
    """Demonstrate CIF file support."""
    print("=" * 60)
    print("CIF File Analysis")
    print("=" * 60)
    
    analyzer = TCRpMHCAnalyzer(verbose=True)
    
    # CIF files are automatically detected and handled
    cif_file = "test_data/structure.cif"  # Replace with your file
    
    try:
        results, chain_map = analyzer.analyze_pdb(cif_file)
        
        print(f"\nCIF file analyzed successfully")
        
        # Chain ID mapping is provided for CIF files
        if chain_map:
            print("\nChain ID mapping (internal -> author):")
            for internal_id, author_id in chain_map.items():
                print(f"  {internal_id} -> {author_id}")
        
        print(f"\nFound {len(results)} complex(es)")
    
    except FileNotFoundError:
        print(f"File {cif_file} not found. Please provide a valid CIF file.")
    except Exception as e:
        print(f"Error: {e}")


def process_results_programmatically():
    """Example of processing results for downstream analysis."""
    print("=" * 60)
    print("Programmatic Result Processing")
    print("=" * 60)
    
    analyzer = TCRpMHCAnalyzer(verbose=False)
    pdb_file = "test_data/1oga.pdb"  # Replace with your file
    
    try:
        results, chain_map = analyzer.analyze_pdb(pdb_file)
        
        # Extract specific information
        mhc_i_complexes = []
        mhc_ii_complexes = []
        unpaired_tcrs = []
        
        for complex_name, complex_data in results.items():
            for pair in complex_data['pairs']:
                if pair['type'] == 'MHC-I':
                    mhc_i_complexes.append(pair)
                elif pair['type'] == 'MHC-II':
                    mhc_ii_complexes.append(pair)
                elif pair['type'] == 'UNPAIRED_TCR':
                    unpaired_tcrs.append(pair)
        
        # Generate summary report
        print("\nSummary Report:")
        print(f"  MHC-I complexes: {len(mhc_i_complexes)}")
        print(f"  MHC-II complexes: {len(mhc_ii_complexes)}")
        print(f"  Unpaired TCRs: {len(unpaired_tcrs)}")
        
        # Extract peptide sequences (example)
        peptides = []
        for complex_name, complex_data in results.items():
            for chain_id, chain_type in complex_data['chains'].items():
                if chain_type == 'PEPTIDE':
                    peptides.append(chain_id)
        
        if peptides:
            print(f"\nPeptide chains found: {', '.join(peptides)}")
    
    except Exception as e:
        print(f"Error: {e}")


if __name__ == "__main__":
    # Run examples
    print("\n" + "=" * 60)
    print("TCRpMHC Analyzer - Basic Usage Examples")
    print("=" * 60 + "\n")
    
    # Note: Replace file paths with actual PDB/CIF files
    print("Note: These examples use placeholder file paths.")
    print("Replace 'test_data/*.pdb' with your actual PDB/CIF files.\n")
    
    # Uncomment to run examples:
    # analyze_single_file()
    # analyze_with_custom_parameters()
    # analyze_cif_file()
    # process_results_programmatically()
    
    print("\nTo run examples, uncomment the function calls at the end of the script.")