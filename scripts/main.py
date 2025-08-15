#!/usr/bin/env python3
"""
pHLA-TCR Splicer - Main Entry Point
===================================

Unified command-line interface for analyzing TCR-pMHC protein complex structures.
Automatically handles both MHC-I (with B2M) and MHC-II (α/β heterodimer) complexes
using sequence alignment, spatial clustering, and pattern recognition.

Usage:
    python main.py <pdb_file>                      # Analyze single file
    python main.py *.pdb --batch-summary           # Batch analysis
    python main.py <pdb_file> --verbose            # Verbose output
    python main.py <pdb_file> --output results.txt # Save to file

Advanced options:
    python main.py <pdb_file> --eps 45.0           # Tighter clustering
    python main.py <pdb_file> --align-score 25     # Stricter TCR alignment

Examples:
    python main.py 1oga.pdb
    python main.py 4z7u.cif --verbose
    python main.py data/*.pdb --batch-summary --output-dir results/

Author: skblnw
"""

import argparse
import sys
import os
import csv
from pathlib import Path
from typing import Dict, List
import logging

# Import our unified analyzer
try:
    from phlatcr_splicer import TCRpMHCAnalyzer
except ImportError:
    # Try adding current directory to path for development
    sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
    try:
        from phlatcr_splicer import TCRpMHCAnalyzer
    except ImportError:
        print("Error: Could not import phlatcr_splicer package.")
        print("Please run from the project root directory or install with: pip install -e .")
        sys.exit(1)

# Configure logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s'
)
logger = logging.getLogger(__name__)


def format_chain_type(chain_type: str) -> str:
    """Format chain type for display."""
    type_map = {
        'TCR_ALPHA': 'TCR-α',
        'TCR_BETA': 'TCR-β',
        'MHC_I_ALPHA': 'MHC-I Heavy Chain',
        'B2M': 'β2-microglobulin',
        'MHC_II_ALPHA': 'MHC-II α',
        'MHC_II_BETA': 'MHC-II β',
        'PEPTIDE': 'Peptide',
        'UNKNOWN': 'Unknown'
    }
    return type_map.get(chain_type, chain_type)


def format_complex(complex_data: Dict) -> str:
    """Format a complex for display."""
    complex_type = complex_data.get('type', 'Unknown')
    
    if complex_type == 'MHC-I':
        return (f"MHC-I Complex: "
                f"Heavy[{complex_data.get('mhc_alpha', '-')}] + "
                f"B2M[{complex_data.get('b2m', '-')}] + "
                f"Peptide[{complex_data.get('peptide', '-')}]")
    elif complex_type == 'MHC-II':
        return (f"MHC-II Complex: "
                f"α[{complex_data.get('mhc_alpha', '-')}] + "
                f"β[{complex_data.get('mhc_beta', '-')}] + "
                f"Peptide[{complex_data.get('peptide', '-')}]")
    elif complex_type == 'UNPAIRED_TCR':
        return (f"Unpaired TCR: "
                f"α[{complex_data.get('tcr_alpha', '-')}] + "
                f"β[{complex_data.get('tcr_beta', '-')}]")
    else:
        return f"Unknown Complex Type: {complex_type}"


def run_analysis(pdb_file: str, analyzer: TCRpMHCAnalyzer, verbose: bool = False) -> Dict:
    """
    Run analysis on a single PDB file.
    
    Returns:
        Dictionary with analysis results
    """
    if verbose:
        print(f"Analyzing {os.path.basename(pdb_file)}...")
        print("=" * 60)
    
    try:
        results, chain_id_map = analyzer.analyze_pdb(pdb_file)
        
        if not results:
            print(f"Warning: No complexes found in {os.path.basename(pdb_file)}")
            return {}
        
        # Process results
        total_chains = 0
        total_complexes = 0
        mhc_i_complexes = 0
        mhc_ii_complexes = 0
        unpaired_tcrs = 0
        
        for complex_name, complex_data in results.items():
            chains = complex_data.get('chains', {})
            pairs = complex_data.get('pairs', [])
            
            total_chains += len(chains)
            total_complexes += len(pairs)
            
            for pair in pairs:
                if pair.get('type') == 'MHC-I':
                    mhc_i_complexes += 1
                elif pair.get('type') == 'MHC-II':
                    mhc_ii_complexes += 1
                elif pair.get('type') == 'UNPAIRED_TCR':
                    unpaired_tcrs += 1
        
        # Display results
        print(f"\nAnalysis Results for {os.path.basename(pdb_file)}:")
        print("-" * 50)
        
        for complex_name, complex_data in sorted(results.items()):
            print(f"\n{complex_name.replace('_', ' ').title()}:")
            
            # Show chain assignments
            chains = complex_data.get('chains', {})
            if chains:
                print("  Chain Assignments:")
                for chain_id, chain_type in sorted(chains.items()):
                    # Map to author chain ID if available
                    display_id = chain_id_map.get(chain_id, chain_id) if chain_id_map else chain_id
                    print(f"    Chain {display_id}: {format_chain_type(chain_type)}")
            
            # Show paired complexes
            pairs = complex_data.get('pairs', [])
            if pairs:
                print("  Identified Complexes:")
                for i, pair in enumerate(pairs, 1):
                    # Handle TCR chains
                    tcr_info = []
                    if pair.get('tcr_alpha'):
                        tcr_info.append(f"TCR-α[{chain_id_map.get(pair['tcr_alpha'], pair['tcr_alpha']) if chain_id_map else pair['tcr_alpha']}]")
                    if pair.get('tcr_beta'):
                        tcr_info.append(f"TCR-β[{chain_id_map.get(pair['tcr_beta'], pair['tcr_beta']) if chain_id_map else pair['tcr_beta']}]")
                    
                    tcr_str = " + ".join(tcr_info) if tcr_info else "No TCR"
                    
                    # Handle pMHC complex
                    if pair.get('type') == 'MHC-I':
                        mhc_str = (f"MHC-I[{chain_id_map.get(pair.get('mhc_alpha', '-'), pair.get('mhc_alpha', '-')) if chain_id_map else pair.get('mhc_alpha', '-')}] + "
                                  f"B2M[{chain_id_map.get(pair.get('b2m', '-'), pair.get('b2m', '-')) if chain_id_map else pair.get('b2m', '-')}] + "
                                  f"Pep[{chain_id_map.get(pair.get('peptide', '-'), pair.get('peptide', '-')) if chain_id_map else pair.get('peptide', '-')}]")
                    elif pair.get('type') == 'MHC-II':
                        mhc_str = (f"MHC-II-α[{chain_id_map.get(pair.get('mhc_alpha', '-'), pair.get('mhc_alpha', '-')) if chain_id_map else pair.get('mhc_alpha', '-')}] + "
                                  f"MHC-II-β[{chain_id_map.get(pair.get('mhc_beta', '-'), pair.get('mhc_beta', '-')) if chain_id_map else pair.get('mhc_beta', '-')}] + "
                                  f"Pep[{chain_id_map.get(pair.get('peptide', '-'), pair.get('peptide', '-')) if chain_id_map else pair.get('peptide', '-')}]")
                    elif pair.get('type') == 'UNPAIRED_TCR':
                        mhc_str = "No pMHC (Unpaired TCR)"
                    else:
                        mhc_str = "Unknown complex"
                    
                    print(f"    Complex #{i}: {tcr_str} :: {mhc_str}")
        
        # Summary
        print(f"\nSummary:")
        print(f"  Total chains identified: {total_chains}")
        print(f"  Total complexes found: {total_complexes}")
        if mhc_i_complexes:
            print(f"  MHC-I complexes: {mhc_i_complexes}")
        if mhc_ii_complexes:
            print(f"  MHC-II complexes: {mhc_ii_complexes}")
        if unpaired_tcrs:
            print(f"  Unpaired TCRs: {unpaired_tcrs}")
        
        return results
        
    except Exception as e:
        logger.error(f"Error analyzing {pdb_file}: {e}")
        if verbose:
            import traceback
            traceback.print_exc()
        return {}


def write_summary_csv(all_results: List[Dict], output_file: str):
    """Write a summary CSV file for batch processing."""
    if not all_results:
        logger.info("No results to write to CSV")
        return
    
    headers = [
        'PDB File', 'Complex ID', 'Complex Type',
        'TCR_Alpha', 'TCR_Beta',
        'MHC_Alpha', 'MHC_Beta', 'B2M',
        'Peptide'
    ]
    
    rows = []
    for result in all_results:
        pdb_file = result['pdb_file']
        results = result['results']
        
        for complex_name, complex_data in results.items():
            pairs = complex_data.get('pairs', [])
            for i, pair in enumerate(pairs):
                row = {
                    'PDB File': os.path.basename(pdb_file),
                    'Complex ID': f"{complex_name}_{i+1}",
                    'Complex Type': pair.get('type', 'Unknown'),
                    'TCR_Alpha': pair.get('tcr_alpha', ''),
                    'TCR_Beta': pair.get('tcr_beta', ''),
                    'MHC_Alpha': pair.get('mhc_alpha', ''),
                    'MHC_Beta': pair.get('mhc_beta', ''),
                    'B2M': pair.get('b2m', ''),
                    'Peptide': pair.get('peptide', '')
                }
                rows.append(row)
    
    with open(output_file, 'w', newline='', encoding='utf-8') as f:
        writer = csv.DictWriter(f, fieldnames=headers)
        writer.writeheader()
        writer.writerows(rows)
    
    logger.info(f"Summary CSV saved to {output_file}")


def main():
    parser = argparse.ArgumentParser(
        description="Analyze TCR-pMHC structures in PDB/CIF files",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=__doc__
    )
    
    # Required arguments
    parser.add_argument('input', nargs='+', 
                       help='Input PDB/CIF file(s) or directory')
    
    # Output options
    parser.add_argument('-o', '--output', 
                       help='Output file for single analysis')
    parser.add_argument('--output-dir', default='output',
                       help='Output directory for batch analysis (default: output)')
    parser.add_argument('--batch-summary', action='store_true',
                       help='Generate summary CSV for batch processing')
    
    # Analysis options
    parser.add_argument('-v', '--verbose', action='store_true',
                       help='Enable verbose output')
    
    # Tuning parameters
    parser.add_argument('--eps', type=float, default=60.0,
                       help='DBSCAN clustering distance in Å (default: 60.0)')
    parser.add_argument('--align-score', type=int, default=20,
                       help='Min alignment score for TCR classification (default: 20)')
    parser.add_argument('--align-ratio', type=float, default=1.5,
                       help='Score ratio for TCR α/β discrimination (default: 1.5)')
    
    # Distance thresholds
    parser.add_argument('--tcr-pair-dist', type=float, default=50.0,
                       help='Max distance for TCR α/β pairing (default: 50.0)')
    parser.add_argument('--mhc1-pair-dist', type=float, default=50.0,
                       help='Max distance for MHC-I/B2M pairing (default: 50.0)')
    parser.add_argument('--mhc2-pair-dist', type=float, default=50.0,
                       help='Max distance for MHC-II α/β pairing (default: 50.0)')
    parser.add_argument('--pep-mhc1-dist', type=float, default=40.0,
                       help='Max distance for Peptide/MHC-I pairing (default: 40.0)')
    parser.add_argument('--pep-mhc2-dist', type=float, default=60.0,
                       help='Max distance for Peptide/MHC-II pairing (default: 60.0)')
    parser.add_argument('--tcr-pmhc-dist', type=float, default=150.0,
                       help='Max distance for TCR/pMHC pairing (default: 150.0)')
    parser.add_argument('--tcr-pmhc-dist-unpaired', type=float, default=120.0,
                       help='Max distance for single-chain TCR/pMHC (default: 120.0)')
    parser.add_argument('--reclassify-dist', type=float, default=50.0,
                       help='Max distance for reclassifying UNKNOWN chains (default: 50.0)')
    
    args = parser.parse_args()
    
    if args.verbose:
        logging.getLogger().setLevel(logging.DEBUG)
    
    # Initialize analyzer with parameters
    analyzer = TCRpMHCAnalyzer(
        verbose=args.verbose,
        eps=args.eps,
        align_score=args.align_score,
        align_ratio=args.align_ratio,
        tcr_pair_dist=args.tcr_pair_dist,
        mhc1_pair_dist=args.mhc1_pair_dist,
        mhc2_pair_dist=args.mhc2_pair_dist,
        pep_mhc1_dist=args.pep_mhc1_dist,
        pep_mhc2_dist=args.pep_mhc2_dist,
        tcr_pmhc_dist=args.tcr_pmhc_dist,
        tcr_pmhc_dist_unpaired=args.tcr_pmhc_dist_unpaired,
        reclassify_dist=args.reclassify_dist
    )
    
    # Collect input files
    input_files = []
    for input_path in args.input:
        if os.path.isdir(input_path):
            # Process directory
            path = Path(input_path)
            input_files.extend(path.glob('*.pdb'))
            input_files.extend(path.glob('*.cif'))
        else:
            # Process file pattern or single file
            from glob import glob
            matching = glob(input_path)
            input_files.extend([Path(f) for f in matching])
    
    if not input_files:
        print(f"Error: No PDB/CIF files found in: {args.input}")
        sys.exit(1)
    
    print(f"Found {len(input_files)} file(s) to analyze")
    
    # Process files
    all_results = []
    for pdb_file in input_files:
        results = run_analysis(str(pdb_file), analyzer, args.verbose)
        if results:
            all_results.append({
                'pdb_file': str(pdb_file),
                'results': results
            })
        
        # Save individual output if requested
        if args.output and len(input_files) == 1:
            with open(args.output, 'w') as f:
                f.write(f"Analysis Results for {pdb_file.name}\n")
                f.write("=" * 60 + "\n\n")
                for complex_name, complex_data in results.items():
                    f.write(f"{complex_name}:\n")
                    f.write("  Chains:\n")
                    for chain_id, chain_type in complex_data.get('chains', {}).items():
                        f.write(f"    {chain_id}: {chain_type}\n")
                    f.write("  Complexes:\n")
                    for i, pair in enumerate(complex_data.get('pairs', []), 1):
                        f.write(f"    Complex #{i}: {pair}\n")
                    f.write("\n")
            print(f"Results saved to {args.output}")
    
    # Generate batch summary if requested
    if args.batch_summary and len(input_files) > 1:
        os.makedirs(args.output_dir, exist_ok=True)
        summary_file = os.path.join(args.output_dir, 'chain_summary.csv')
        write_summary_csv(all_results, summary_file)
        print(f"Batch summary saved to {summary_file}")
    
    print(f"\nAnalysis complete! Processed {len(input_files)} file(s)")


if __name__ == "__main__":
    main()