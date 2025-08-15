#!/usr/bin/env python3
"""
Batch processing examples for the Unified TCR-pMHC Analyzer

This script demonstrates how to analyze multiple PDB/CIF files
and generate comprehensive summaries.

Author: skblnw
"""

import os
import csv
from pathlib import Path
from typing import List, Dict
from phlatcr_splicer import TCRpMHCAnalyzer


def batch_analyze_directory(directory: str, file_pattern: str = "*.pdb"):
    """
    Analyze all PDB/CIF files in a directory.
    
    Args:
        directory: Path to directory containing structure files
        file_pattern: Glob pattern for files (default: "*.pdb")
    """
    print("=" * 60)
    print(f"Batch Analysis of {directory}")
    print("=" * 60)
    
    # Initialize analyzer
    analyzer = TCRpMHCAnalyzer(verbose=False)
    
    # Find all matching files
    dir_path = Path(directory)
    pdb_files = list(dir_path.glob(file_pattern))
    cif_files = list(dir_path.glob("*.cif")) if file_pattern == "*.pdb" else []
    all_files = pdb_files + cif_files
    
    if not all_files:
        print(f"No files found matching pattern in {directory}")
        return
    
    print(f"Found {len(all_files)} file(s) to analyze\n")
    
    # Store results for summary
    all_results = []
    successful = 0
    failed = 0
    
    # Process each file
    for i, file_path in enumerate(all_files, 1):
        print(f"Processing {i}/{len(all_files)}: {file_path.name}")
        
        try:
            results, chain_map = analyzer.analyze_pdb(str(file_path))
            
            # Count statistics
            total_chains = sum(len(c['chains']) for c in results.values())
            total_complexes = sum(len(c['pairs']) for c in results.values())
            
            print(f"  ✓ Found {total_chains} chains, {total_complexes} complexes")
            
            all_results.append({
                'file': file_path.name,
                'results': results,
                'chain_map': chain_map,
                'status': 'success'
            })
            successful += 1
            
        except Exception as e:
            print(f"  ✗ Error: {e}")
            all_results.append({
                'file': file_path.name,
                'error': str(e),
                'status': 'failed'
            })
            failed += 1
    
    # Print summary
    print("\n" + "=" * 60)
    print("Batch Analysis Summary")
    print("=" * 60)
    print(f"Total files processed: {len(all_files)}")
    print(f"Successful: {successful}")
    print(f"Failed: {failed}")
    
    return all_results


def generate_summary_csv(results: List[Dict], output_file: str):
    """
    Generate a CSV summary of batch analysis results.
    
    Args:
        results: List of analysis results from batch_analyze_directory
        output_file: Path to output CSV file
    """
    print(f"\nGenerating CSV summary: {output_file}")
    
    # Prepare CSV headers
    headers = [
        'File', 'Complex_ID', 'Complex_Type',
        'TCR_Alpha', 'TCR_Beta',
        'MHC_Alpha', 'MHC_Beta', 'B2M',
        'Peptide', 'Status'
    ]
    
    rows = []
    
    for result in results:
        if result['status'] == 'failed':
            # Add error row
            rows.append({
                'File': result['file'],
                'Status': f"Error: {result.get('error', 'Unknown')}"
            })
            continue
        
        # Process successful results
        for complex_name, complex_data in result['results'].items():
            for i, pair in enumerate(complex_data['pairs']):
                row = {
                    'File': result['file'],
                    'Complex_ID': f"{complex_name}_{i+1}",
                    'Complex_Type': pair.get('type', 'Unknown'),
                    'TCR_Alpha': pair.get('tcr_alpha', ''),
                    'TCR_Beta': pair.get('tcr_beta', ''),
                    'MHC_Alpha': pair.get('mhc_alpha', ''),
                    'MHC_Beta': pair.get('mhc_beta', ''),
                    'B2M': pair.get('b2m', ''),
                    'Peptide': pair.get('peptide', ''),
                    'Status': 'Success'
                }
                rows.append(row)
    
    # Write CSV
    with open(output_file, 'w', newline='', encoding='utf-8') as f:
        writer = csv.DictWriter(f, fieldnames=headers)
        writer.writeheader()
        writer.writerows(rows)
    
    print(f"CSV summary saved with {len(rows)} entries")


def analyze_with_filtering(directory: str, min_chains: int = 3):
    """
    Analyze files and filter results based on criteria.
    
    Args:
        directory: Directory containing PDB files
        min_chains: Minimum number of chains required
    """
    print("=" * 60)
    print("Filtered Batch Analysis")
    print("=" * 60)
    
    analyzer = TCRpMHCAnalyzer(verbose=False)
    
    # Find files
    dir_path = Path(directory)
    all_files = list(dir_path.glob("*.pdb")) + list(dir_path.glob("*.cif"))
    
    print(f"Analyzing files with >= {min_chains} chains\n")
    
    filtered_results = []
    
    for file_path in all_files:
        try:
            results, chain_map = analyzer.analyze_pdb(str(file_path))
            
            # Apply filter
            total_chains = sum(len(c['chains']) for c in results.values())
            
            if total_chains >= min_chains:
                print(f"✓ {file_path.name}: {total_chains} chains (included)")
                filtered_results.append({
                    'file': file_path.name,
                    'chains': total_chains,
                    'results': results
                })
            else:
                print(f"✗ {file_path.name}: {total_chains} chains (excluded)")
        
        except Exception as e:
            print(f"⚠ {file_path.name}: Error - {e}")
    
    print(f"\n{len(filtered_results)} files passed the filter")
    return filtered_results


def parallel_batch_processing(file_list: List[str], num_workers: int = 4):
    """
    Process multiple files in parallel for better performance.
    
    Args:
        file_list: List of PDB/CIF file paths
        num_workers: Number of parallel workers
    """
    print("=" * 60)
    print(f"Parallel Batch Processing ({num_workers} workers)")
    print("=" * 60)
    
    from concurrent.futures import ProcessPoolExecutor, as_completed
    
    def analyze_file(file_path):
        """Worker function to analyze a single file."""
        analyzer = TCRpMHCAnalyzer(verbose=False)
        try:
            results, chain_map = analyzer.analyze_pdb(file_path)
            return {
                'file': os.path.basename(file_path),
                'results': results,
                'chain_map': chain_map,
                'status': 'success'
            }
        except Exception as e:
            return {
                'file': os.path.basename(file_path),
                'error': str(e),
                'status': 'failed'
            }
    
    print(f"Processing {len(file_list)} files...\n")
    
    all_results = []
    
    with ProcessPoolExecutor(max_workers=num_workers) as executor:
        # Submit all tasks
        future_to_file = {
            executor.submit(analyze_file, f): f 
            for f in file_list
        }
        
        # Process completed tasks
        for future in as_completed(future_to_file):
            file_path = future_to_file[future]
            try:
                result = future.result()
                status = "✓" if result['status'] == 'success' else "✗"
                print(f"{status} {os.path.basename(file_path)}")
                all_results.append(result)
            except Exception as e:
                print(f"✗ {os.path.basename(file_path)}: {e}")
    
    # Summary
    successful = sum(1 for r in all_results if r['status'] == 'success')
    failed = len(all_results) - successful
    
    print(f"\nCompleted: {successful} successful, {failed} failed")
    return all_results


def generate_statistics_report(results: List[Dict]):
    """
    Generate detailed statistics from batch analysis results.
    
    Args:
        results: List of analysis results
    """
    print("=" * 60)
    print("Statistical Analysis Report")
    print("=" * 60)
    
    # Collect statistics
    stats = {
        'total_files': len(results),
        'successful': 0,
        'failed': 0,
        'total_chains': 0,
        'total_complexes': 0,
        'mhc_i_complexes': 0,
        'mhc_ii_complexes': 0,
        'unpaired_tcrs': 0,
        'chain_types': {},
        'files_with_multiple_complexes': 0
    }
    
    for result in results:
        if result['status'] == 'failed':
            stats['failed'] += 1
            continue
        
        stats['successful'] += 1
        
        # Count complexes in this file
        file_complexes = 0
        
        for complex_data in result['results'].values():
            # Count chains
            for chain_type in complex_data['chains'].values():
                stats['chain_types'][chain_type] = stats['chain_types'].get(chain_type, 0) + 1
                stats['total_chains'] += 1
            
            # Count complex types
            for pair in complex_data['pairs']:
                file_complexes += 1
                stats['total_complexes'] += 1
                
                if pair['type'] == 'MHC-I':
                    stats['mhc_i_complexes'] += 1
                elif pair['type'] == 'MHC-II':
                    stats['mhc_ii_complexes'] += 1
                elif pair['type'] == 'UNPAIRED_TCR':
                    stats['unpaired_tcrs'] += 1
        
        if file_complexes > 1:
            stats['files_with_multiple_complexes'] += 1
    
    # Print report
    print(f"\nFiles Analyzed: {stats['total_files']}")
    print(f"  Successful: {stats['successful']}")
    print(f"  Failed: {stats['failed']}")
    
    print(f"\nChains Identified: {stats['total_chains']}")
    print("  Chain Type Distribution:")
    for chain_type, count in sorted(stats['chain_types'].items()):
        print(f"    {chain_type}: {count}")
    
    print(f"\nComplexes Found: {stats['total_complexes']}")
    print(f"  MHC-I complexes: {stats['mhc_i_complexes']}")
    print(f"  MHC-II complexes: {stats['mhc_ii_complexes']}")
    print(f"  Unpaired TCRs: {stats['unpaired_tcrs']}")
    
    print(f"\nFiles with multiple complexes: {stats['files_with_multiple_complexes']}")
    
    return stats


if __name__ == "__main__":
    print("\n" + "=" * 60)
    print("TCRpMHC Analyzer - Batch Processing Examples")
    print("=" * 60 + "\n")
    
    # Example usage (replace with actual paths)
    test_directory = "test_data"
    
    # Example 1: Basic batch analysis
    # results = batch_analyze_directory(test_directory)
    
    # Example 2: Generate CSV summary
    # if results:
    #     generate_summary_csv(results, "batch_summary.csv")
    
    # Example 3: Filtered analysis
    # filtered = analyze_with_filtering(test_directory, min_chains=5)
    
    # Example 4: Parallel processing (for many files)
    # file_list = ["file1.pdb", "file2.pdb", "file3.pdb"]
    # parallel_results = parallel_batch_processing(file_list)
    
    # Example 5: Generate statistics
    # if results:
    #     stats = generate_statistics_report(results)
    
    print("Note: Uncomment the example function calls to run them.")
    print("Replace 'test_data' with your actual directory path.")