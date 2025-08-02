#!/usr/bin/env python3
"""
pHLA-TCR Splicer - Main Entry Point
===================================

A unified command-line interface for analyzing protein complex structures.
Supports both MHC-I (pHLA-TCR) and MHC-II (pMHC-II-TCR) complexes.

Usage:
    python main.py --type mhc-i <pdb_file>     # For MHC-I complexes
    python main.py --type mhc-ii <pdb_file>    # For MHC-II complexes
    python main.py --auto <pdb_file>           # Auto-detect complex type

Examples:
    python main.py --type mhc-i 1oga.pdb
    python main.py --type mhc-ii 4z7u.pdb
    python main.py --auto complex.pdb --verbose --output results.txt
"""

import argparse
import sys
import os
from pathlib import Path

# Import our analyzers
try:
    from phlatcr_splicer import pMHCITCRAnalyzer, pMHCIITCRAnalyzer
except ImportError:
    # Try adding current directory to path for development
    sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
    try:
        from phlatcr_splicer import pMHCITCRAnalyzer, pMHCIITCRAnalyzer
    except ImportError:
        print("Error: Could not import phlatcr_splicer package.")
        print("Please run from the project root directory or install with: pip install -e .")
        sys.exit(1)


def detect_complex_type(pdb_file: str, verbose: bool = False) -> str:
    """
    Auto-detect whether a PDB file contains MHC-I or MHC-II complexes.
    
    Strategy:
    1. Look for β2-microglobulin (indicates MHC-I)
    2. Count peptide lengths (MHC-I: 8-11, MHC-II: 12-25)
    3. Analyze chain count and composition
    """
    if verbose:
        print("🔍 Auto-detecting complex type...")
    
    # Try both analyzers quickly to see which gives better results
    try:
        # Test MHC-I analyzer
        mhc_i_analyzer = pMHCITCRAnalyzer(verbose=False)
        mhc_i_result = mhc_i_analyzer.analyze_pdb(pdb_file)
        
        # Test MHC-II analyzer  
        mhc_ii_analyzer = pMHCIITCRAnalyzer(verbose=False)
        mhc_ii_result = mhc_ii_analyzer.analyze_pdb(pdb_file)
        
        # Count successful identifications and unknown chains
        mhc_i_unknown = sum(1 for chain_type in mhc_i_result.values() if chain_type == 'unknown')
        mhc_ii_unknown = sum(1 for chain_type in mhc_ii_result.values() if chain_type == 'unknown')
        
        # Check for b2m (strong indicator of MHC-I)
        has_b2m = any('b2m' in chain_type for chain_type in mhc_i_result.values())
        
        # Check for MHC-II specific patterns
        has_mhc_ii_alpha = any('mhc_ii_alpha' in chain_type for chain_type in mhc_ii_result.values())
        has_mhc_ii_beta = any('mhc_ii_beta' in chain_type for chain_type in mhc_ii_result.values())
        
        # Check for peptide length indicators
        has_long_peptides = any('peptide' in chain_type for chain_type in mhc_ii_result.values())
        has_short_peptides = any('peptide' in chain_type for chain_type in mhc_i_result.values())
        
        # Decision logic with improved scoring
        if has_b2m:
            detected_type = "mhc-i"
            reason = "β2-microglobulin detected"
        elif has_mhc_ii_alpha and has_mhc_ii_beta:
            detected_type = "mhc-ii"
            reason = "MHC-II α/β heterodimer detected"
        elif mhc_ii_unknown == 0 and mhc_i_unknown > 0:
            detected_type = "mhc-ii"
            reason = "MHC-II analysis perfect, MHC-I has unknown chains"
        elif mhc_i_unknown == 0 and mhc_ii_unknown > 0:
            detected_type = "mhc-i"
            reason = "MHC-I analysis perfect, MHC-II has unknown chains"
        elif mhc_ii_unknown < mhc_i_unknown:
            detected_type = "mhc-ii"
            reason = f"fewer unknown chains with MHC-II analysis ({mhc_ii_unknown} vs {mhc_i_unknown})"
        elif mhc_i_unknown < mhc_ii_unknown:
            detected_type = "mhc-i"
            reason = f"fewer unknown chains with MHC-I analysis ({mhc_i_unknown} vs {mhc_ii_unknown})"
        elif len(mhc_ii_result) >= 8:  # Large number of chains suggests MHC-II multi-complex
            detected_type = "mhc-ii"
            reason = "large number of chains suggests MHC-II multi-complex"
        else:
            detected_type = "mhc-i"  # Default fallback
            reason = "defaulting to MHC-I (both analyses similar)"
        
        if verbose:
            print(f"✅ Detected: {detected_type.upper()} complex ({reason})")
            print(f"   MHC-I analysis: {len(mhc_i_result) - mhc_i_unknown}/{len(mhc_i_result)} chains identified")
            print(f"   MHC-II analysis: {len(mhc_ii_result) - mhc_ii_unknown}/{len(mhc_ii_result)} chains identified")
        
        return detected_type
        
    except Exception as e:
        if verbose:
            print(f"⚠️  Auto-detection failed: {e}")
            print("🔄 Defaulting to MHC-I analyzer")
        return "mhc-i"


def run_analysis(analyzer_type: str, pdb_file: str, verbose: bool = False, output_file: str = None) -> dict:
    """Run the specified analyzer on the PDB file."""
    
    # Initialize the appropriate analyzer
    if analyzer_type == "mhc-i":
        analyzer = pMHCITCRAnalyzer(verbose=verbose)
        analyzer_name = "MHC-I (pHLA-TCR)"
    elif analyzer_type == "mhc-ii":
        analyzer = pMHCIITCRAnalyzer(verbose=verbose)
        analyzer_name = "MHC-II (pMHC-II-TCR)"
    else:
        raise ValueError(f"Unknown analyzer type: {analyzer_type}")
    
    print(f"🧬 Running {analyzer_name} Analysis")
    print("=" * 50)
    
    # Run analysis
    try:
        results = analyzer.analyze_pdb(pdb_file)
        
        # Display results
        print(f"\n📊 Analysis Results for {os.path.basename(pdb_file)}:")
        print("-" * 40)
        
        for chain_id in sorted(results.keys()):
            chain_type = results[chain_id]
            print(f"  Chain {chain_id}: {chain_type}")
        
        # Summary statistics
        total_chains = len(results)
        unknown_chains = sum(1 for chain_type in results.values() if chain_type == 'unknown')
        identified_chains = total_chains - unknown_chains
        
        print(f"\n✅ Summary:")
        print(f"  Total chains: {total_chains}")
        print(f"  Identified: {identified_chains}")
        print(f"  Unknown: {unknown_chains}")
        
        if unknown_chains == 0:
            print("🎯 Perfect! All chains identified successfully.")
        elif unknown_chains <= 2:
            print("👍 Good analysis with minimal unknown chains.")
        else:
            print("⚠️  Some chains could not be identified. Consider trying the other analyzer.")
        
        # Save to file if requested
        if output_file:
            save_results_to_file(results, pdb_file, analyzer_name, output_file, verbose)
        
        return results
        
    except FileNotFoundError:
        print(f"❌ Error: PDB file '{pdb_file}' not found.")
        sys.exit(1)
    except Exception as e:
        print(f"❌ Error during analysis: {e}")
        sys.exit(1)


def save_results_to_file(results: dict, pdb_file: str, analyzer_name: str, output_file: str, verbose: bool = False):
    """Save analysis results to a file."""
    try:
        with open(output_file, 'w') as f:
            f.write(f"pHLA-TCR Splicer Analysis Report\n")
            f.write(f"================================\n\n")
            f.write(f"Analyzer: {analyzer_name}\n")
            f.write(f"PDB File: {pdb_file}\n")
            f.write(f"Timestamp: {__import__('datetime').datetime.now().strftime('%Y-%m-%d %H:%M:%S')}\n\n")
            
            f.write("Chain Assignments:\n")
            f.write("-" * 18 + "\n")
            for chain_id in sorted(results.keys()):
                chain_type = results[chain_id]
                f.write(f"Chain {chain_id}: {chain_type}\n")
            
            # Summary
            total_chains = len(results)
            unknown_chains = sum(1 for chain_type in results.values() if chain_type == 'unknown')
            f.write(f"\nSummary:\n")
            f.write(f"Total chains: {total_chains}\n")
            f.write(f"Identified: {total_chains - unknown_chains}\n")
            f.write(f"Unknown: {unknown_chains}\n")
        
        if verbose:
            print(f"💾 Results saved to: {output_file}")
            
    except Exception as e:
        print(f"⚠️  Warning: Could not save results to {output_file}: {e}")


def main():
    """Main entry point for the pHLA-TCR Splicer."""
    parser = argparse.ArgumentParser(
        description="pHLA-TCR Splicer: Analyze protein complex structures",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  %(prog)s --type mhc-i 1oga.pdb
  %(prog)s --type mhc-ii 4z7u.pdb  
  %(prog)s --auto complex.pdb --verbose
  %(prog)s --type mhc-i input.pdb --output results.txt

Analyzer Types:
  mhc-i   : MHC Class I complexes (pHLA-TCR)
            - MHC heavy chain, β2-microglobulin, short peptides (8-11), TCR α/β
  
  mhc-ii  : MHC Class II complexes (pMHC-II-TCR)  
            - MHC-II α/β chains, longer peptides (12-25), TCR α/β
  
  auto    : Automatically detect complex type (experimental)
        """
    )
    
    # Analyzer selection (mutually exclusive)
    analyzer_group = parser.add_mutually_exclusive_group(required=True)
    analyzer_group.add_argument(
        '--type', 
        choices=['mhc-i', 'mhc-ii'],
        help='Specify analyzer type'
    )
    analyzer_group.add_argument(
        '--auto', 
        action='store_true',
        help='Auto-detect complex type (experimental)'
    )
    
    # Required PDB file
    parser.add_argument(
        'pdb_file',
        help='Path to PDB file to analyze'
    )
    
    # Optional arguments
    parser.add_argument(
        '--output', '-o',
        help='Output file to save results (optional)'
    )
    parser.add_argument(
        '--verbose', '-v',
        action='store_true',
        help='Enable verbose output'
    )
    parser.add_argument(
        '--version',
        action='version',
        version='pHLA-TCR Splicer v1.0.0'
    )
    
    # Parse arguments
    args = parser.parse_args()
    
    # Validate PDB file exists
    if not os.path.exists(args.pdb_file):
        print(f"❌ Error: PDB file '{args.pdb_file}' does not exist.")
        sys.exit(1)
    
    # Determine analyzer type
    if args.auto:
        analyzer_type = detect_complex_type(args.pdb_file, args.verbose)
    else:
        analyzer_type = args.type
    
    # Run analysis
    try:
        results = run_analysis(analyzer_type, args.pdb_file, args.verbose, args.output)
        
        print("\n🎉 Analysis completed successfully!")
        if args.verbose:
            print(f"📁 Analyzed: {args.pdb_file}")
            print(f"🔬 Used: {analyzer_type.upper()} analyzer")
            if args.output:
                print(f"💾 Saved: {args.output}")
        
    except KeyboardInterrupt:
        print("\n⚠️  Analysis interrupted by user.")
        sys.exit(1)


if __name__ == "__main__":
    main()