"""
Dual Analyzer Example - MHC-I and MHC-II Complex Analysis

This example demonstrates how to use both the MHC-I (pHLA-TCR) and MHC-II (pMHC-II-TCR) 
analyzers from the phlatcr_splicer package.

Author: AI Assistant
"""

from phlatcr_splicer import pMHCITCRAnalyzer, pMHCIITCRAnalyzer


def analyze_mhc_i_complex(pdb_file: str):
    """
    Analyze an MHC-I complex (pHLA-TCR).
    
    Expected components:
    - MHC heavy chain (~270-380 residues)
    - β2-microglobulin (~99 residues)  
    - Short peptide (8-11 residues)
    - TCR α chain (180-230 residues)
    - TCR β chain (230-290 residues)
    """
    print(f"🔬 Analyzing MHC-I Complex: {pdb_file}")
    print("-" * 40)
    
    # Initialize MHC-I analyzer
    mhc_i_analyzer = pMHCITCRAnalyzer(verbose=True)
    
    # Analyze the complex
    result = mhc_i_analyzer.analyze_pdb(pdb_file)
    
    # Display results
    print("\n📊 MHC-I Analysis Results:")
    for chain_id in sorted(result.keys()):
        chain_type = result[chain_id]
        print(f"  Chain {chain_id}: {chain_type}")
    
    # Save detailed report
    mhc_i_analyzer.save_analysis_report(result, f"{pdb_file}_mhc_i_report.txt")
    
    return result


def analyze_mhc_ii_complex(pdb_file: str):
    """
    Analyze an MHC-II complex (pMHC-II-TCR).
    
    Expected components:
    - MHC-II α chain (~180-200 residues)
    - MHC-II β chain (~190-210 residues)
    - Longer peptide (12-25 residues)
    - TCR α chain (180-230 residues)  
    - TCR β chain (230-290 residues)
    """
    print(f"🔬 Analyzing MHC-II Complex: {pdb_file}")
    print("-" * 40)
    
    # Initialize MHC-II analyzer
    mhc_ii_analyzer = pMHCIITCRAnalyzer(verbose=True)
    
    # Analyze the complex
    result = mhc_ii_analyzer.analyze_pdb(pdb_file)
    
    # Display results
    print("\n📊 MHC-II Analysis Results:")
    for chain_id in sorted(result.keys()):
        chain_type = result[chain_id]
        print(f"  Chain {chain_id}: {chain_type}")
    
    # Save detailed report
    mhc_ii_analyzer.save_analysis_report(result, f"{pdb_file}_mhc_ii_report.txt")
    
    return result


def compare_analyzers():
    """
    Compare the capabilities of both analyzers.
    """
    print("🔍 ANALYZER COMPARISON")
    print("=" * 50)
    
    # Initialize both analyzers
    mhc_i = pMHCITCRAnalyzer(verbose=False)
    mhc_ii = pMHCIITCRAnalyzer(verbose=False)
    
    print("\n📋 Chain Type Differences:")
    print("┌─ MHC-I Complexes:")
    print("│  ├─ MHC Heavy Chain (single chain)")
    print("│  ├─ β2-microglobulin (required)")
    print("│  └─ Short peptides (8-11 residues)")
    print("│")
    print("└─ MHC-II Complexes:")
    print("   ├─ MHC-II α + β chains (heterodimer)")
    print("   ├─ No β2-microglobulin")
    print("   └─ Longer peptides (12-25 residues)")
    
    print("\n🧬 Pattern Recognition Examples:")
    print("MHC-I Heavy Chain Signatures:")
    for pattern in mhc_i.mhc_patterns['highly_specific_motifs'][:4]:
        print(f"  • {pattern}")
    
    print("\nMHC-II α Chain Signatures:")
    for pattern in mhc_ii.mhc_ii_patterns['alpha']['highly_specific_motifs'][:4]:
        print(f"  • {pattern}")
    
    print("\nMHC-II β Chain Signatures:")
    for pattern in mhc_ii.mhc_ii_patterns['beta']['highly_specific_motifs'][:4]:
        print(f"  • {pattern}")
    
    print(f"\n📏 Peptide Length Ranges:")
    print(f"  MHC-I:  {mhc_i.peptide_patterns['length_range']} residues")
    print(f"  MHC-II: {mhc_ii.peptide_patterns['length_range']} residues")


def batch_analysis(pdb_files: list, complex_type: str = 'auto'):
    """
    Perform batch analysis on multiple PDB files.
    
    Args:
        pdb_files: List of PDB file paths
        complex_type: 'mhc_i', 'mhc_ii', or 'auto' for automatic detection
    """
    print(f"🚀 BATCH ANALYSIS ({len(pdb_files)} files)")
    print("=" * 50)
    
    results = {}
    
    for pdb_file in pdb_files:
        print(f"\n📁 Processing: {pdb_file}")
        
        if complex_type == 'auto':
            # Try to determine complex type based on file patterns or analysis
            # For demo purposes, we'll analyze with both if requested
            print("🤖 Auto-detection mode: analyzing with both analyzers")
            
            try:
                mhc_i_result = analyze_mhc_i_complex(pdb_file)
                mhc_ii_result = analyze_mhc_ii_complex(pdb_file)
                
                # Simple heuristic: if we find β2m, it's likely MHC-I
                has_b2m = any('b2m' in chain_type for chain_type in mhc_i_result.values())
                has_mhc_ii_chains = any('mhc_ii' in chain_type for chain_type in mhc_ii_result.values())
                
                if has_b2m and not has_mhc_ii_chains:
                    print("  ✅ Determined: MHC-I complex")
                    results[pdb_file] = {'type': 'MHC-I', 'result': mhc_i_result}
                elif has_mhc_ii_chains and not has_b2m:
                    print("  ✅ Determined: MHC-II complex")
                    results[pdb_file] = {'type': 'MHC-II', 'result': mhc_ii_result}
                else:
                    print("  ⚠️  Ambiguous - saved both analyses")
                    results[pdb_file] = {'type': 'Both', 'mhc_i': mhc_i_result, 'mhc_ii': mhc_ii_result}
                    
            except Exception as e:
                print(f"  ❌ Error: {e}")
                results[pdb_file] = {'type': 'Error', 'error': str(e)}
                
        elif complex_type == 'mhc_i':
            try:
                result = analyze_mhc_i_complex(pdb_file)
                results[pdb_file] = {'type': 'MHC-I', 'result': result}
            except Exception as e:
                results[pdb_file] = {'type': 'Error', 'error': str(e)}
                
        elif complex_type == 'mhc_ii':
            try:
                result = analyze_mhc_ii_complex(pdb_file)
                results[pdb_file] = {'type': 'MHC-II', 'result': result}
            except Exception as e:
                results[pdb_file] = {'type': 'Error', 'error': str(e)}
    
    return results


def main():
    """
    Main example function demonstrating both analyzers.
    """
    print("🧬 DUAL ANALYZER DEMONSTRATION")
    print("=" * 60)
    
    # Show analyzer comparison
    compare_analyzers()
    
    print("\n" + "=" * 60)
    print("📝 USAGE EXAMPLES")
    print("=" * 60)
    
    print("\n1. Individual Analysis:")
    print("   # For MHC-I complexes")
    print("   result1 = analyze_mhc_i_complex('1oga.pdb')")
    print("   ")
    print("   # For MHC-II complexes")  
    print("   result2 = analyze_mhc_ii_complex('mhc_ii_complex.pdb')")
    
    print("\n2. Batch Analysis:")
    print("   files = ['complex1.pdb', 'complex2.pdb', 'complex3.pdb']")
    print("   results = batch_analysis(files, complex_type='auto')")
    
    print("\n3. Python API:")
    print("   from phlatcr_splicer import pMHCITCRAnalyzer, pMHCIITCRAnalyzer")
    print("   ")
    print("   # MHC-I analysis")
    print("   mhc_i_analyzer = pMHCITCRAnalyzer(verbose=True)")
    print("   result1 = mhc_i_analyzer.analyze_pdb('complex.pdb')")
    print("   ")
    print("   # MHC-II analysis")
    print("   mhc_ii_analyzer = pMHCIITCRAnalyzer(verbose=True)")
    print("   result2 = mhc_ii_analyzer.analyze_pdb('complex.pdb')")
    
    print("\n4. Command Line:")
    print("   mhc-i-analyze complex.pdb       # MHC-I")
    print("   mhc-ii-analyze complex.pdb      # MHC-II")
    
    print("\n✅ Both analyzers ready for use!")
    print("   📦 Package: phlatcr_splicer")
    print("   🔬 MHC-I: pMHCITCRAnalyzer")
    print("   🔬 MHC-II: pMHCIITCRAnalyzer")


if __name__ == "__main__":
    main()