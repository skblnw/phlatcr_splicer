#!/usr/bin/env python3
"""
Example usage script for pHLA-TCR Complex Structure Analyzer

This script demonstrates how to use the analyzer with various types of input.
"""

import os
import sys
from phlatcr_analyzer import pHLATCRAnalyzer


def example_basic_usage():
    """Basic usage example."""
    print("Example 1: Basic Usage")
    print("-" * 30)
    
    # Initialize the analyzer
    analyzer = pHLATCRAnalyzer(verbose=True)
    
    # Example PDB file (you would replace this with your actual file)
    pdb_file = "example_complex.pdb"
    
    if not os.path.exists(pdb_file):
        print(f"Note: {pdb_file} not found. This is just an example.")
        print("Replace with your actual PDB file path.")
        return
    
    try:
        # Analyze the PDB file
        result = analyzer.analyze_pdb(pdb_file)
        
        # Print results
        print("\nChain assignments:")
        for chain_id, chain_type in result.items():
            print(f"  Chain {chain_id}: {chain_type}")
        
        # Save detailed report
        analyzer.save_analysis_report(result, "analysis_report.txt")
        print("\nDetailed report saved to: analysis_report.txt")
        
    except Exception as e:
        print(f"Error analyzing {pdb_file}: {e}")


def example_batch_processing():
    """Example of processing multiple PDB files."""
    print("\nExample 2: Batch Processing")
    print("-" * 30)
    
    # List of PDB files to process
    pdb_files = [
        "complex1.pdb",
        "complex2.pdb", 
        "complex3.pdb"
    ]
    
    analyzer = pHLATCRAnalyzer(verbose=False)
    
    results = {}
    
    for pdb_file in pdb_files:
        if os.path.exists(pdb_file):
            try:
                result = analyzer.analyze_pdb(pdb_file)
                results[pdb_file] = result
                print(f"✓ Processed {pdb_file}")
            except Exception as e:
                print(f"✗ Error processing {pdb_file}: {e}")
        else:
            print(f"- Skipped {pdb_file} (not found)")
    
    # Summary of all results
    print(f"\nProcessed {len(results)} files successfully")
    
    # Save batch report
    with open("batch_analysis_report.txt", "w") as f:
        f.write("Batch Analysis Report\\n")
        f.write("=" * 50 + "\\n\\n")
        
        for pdb_file, result in results.items():
            f.write(f"File: {pdb_file}\\n")
            for chain_id, chain_type in result.items():
                f.write(f"  Chain {chain_id}: {chain_type}\\n")
            f.write("\\n")


def example_custom_analysis():
    """Example of custom analysis with additional processing."""
    print("\nExample 3: Custom Analysis")
    print("-" * 30)
    
    class CustomAnalyzer(pHLATCRAnalyzer):
        """Extended analyzer with custom functionality."""
        
        def analyze_pdb_detailed(self, pdb_file):
            """Perform detailed analysis with additional metrics."""
            # Get basic analysis
            assignments = self.analyze_pdb(pdb_file)
            
            # Add custom analysis
            structure = self.parser.get_structure('complex', pdb_file)
            chain_info = self._extract_chain_info(structure)
            
            detailed_results = {}
            for chain_id, chain_type in assignments.items():
                info = chain_info[chain_id]
                detailed_results[chain_id] = {
                    'type': chain_type,
                    'length': info['length'],
                    'molecular_weight': info['properties'].get('molecular_weight', 0),
                    'sequence': info['sequence']
                }
            
            return detailed_results
    
    # Use custom analyzer
    analyzer = CustomAnalyzer(verbose=True)
    
    pdb_file = "example_complex.pdb"
    if os.path.exists(pdb_file):
        try:
            detailed_result = analyzer.analyze_pdb_detailed(pdb_file)
            
            print("Detailed analysis results:")
            for chain_id, details in detailed_result.items():
                print(f"\\nChain {chain_id}:")
                print(f"  Type: {details['type']}")
                print(f"  Length: {details['length']} residues")
                print(f"  Molecular Weight: {details['molecular_weight']:.1f} Da")
                print(f"  Sequence: {details['sequence'][:50]}...")
                
        except Exception as e:
            print(f"Error in detailed analysis: {e}")
    else:
        print(f"Example file {pdb_file} not found.")


def create_sample_input():
    """Create sample input files for testing."""
    print("\\nCreating sample input files...")
    
    # Create a minimal sample PDB file
    sample_pdb_content = '''HEADER    IMMUNE SYSTEM/IMMUNE SYSTEM             01-JAN-20   SAMP    
TITLE     SAMPLE PHLA-TCR COMPLEX
ATOM      1  N   MET A   1      20.154  16.772  12.000  1.00 20.00           N  
ATOM      2  CA  MET A   1      19.000  16.500  12.500  1.00 20.00           C  
ATOM      3  C   MET A   1      18.000  15.500  13.000  1.00 20.00           C  
ATOM      4  O   MET A   1      17.500  15.000  12.000  1.00 20.00           O  
ATOM      5  N   ALA A   2      17.800  15.200  14.200  1.00 20.00           N  
ATOM      6  CA  ALA A   2      16.900  14.100  14.700  1.00 20.00           C  
ATOM      7  N   GLY B   1      25.000  20.000  15.000  1.00 20.00           N  
ATOM      8  CA  GLY B   1      24.500  20.500  15.500  1.00 20.00           C  
ATOM      9  C   GLY B   1      23.500  19.500  16.000  1.00 20.00           C  
ATOM     10  O   GLY B   1      23.000  19.000  15.000  1.00 20.00           O  
END
'''
    
    with open("sample_complex.pdb", "w") as f:
        f.write(sample_pdb_content)
    
    print("Created sample_complex.pdb")
    
    # Create usage instructions
    instructions = '''Usage Instructions for pHLA-TCR Analyzer
==========================================

1. Basic command-line usage:
   python phlatcr_analyzer.py your_complex.pdb

2. With verbose output:
   python phlatcr_analyzer.py -v your_complex.pdb

3. Save report to file:
   python phlatcr_analyzer.py -o report.txt your_complex.pdb

4. In Python script:
   from phlatcr_analyzer import pHLATCRAnalyzer
   analyzer = pHLATCRAnalyzer()
   result = analyzer.analyze_pdb("your_complex.pdb")

Expected Chain Types:
- peptide: Short peptide antigen (typically 8-15 residues)
- mhc_heavy: MHC class I heavy chain (270-380 residues)
- b2m: β2-microglobulin (95-105 residues)
- tcr_alpha: TCR α chain (200-250 residues)
- tcr_beta: TCR β chain (240-290 residues)

For more examples, see example_usage.py
'''
    
    with open("USAGE.txt", "w") as f:
        f.write(instructions)
    
    print("Created USAGE.txt")


def main():
    """Run all examples."""
    print("pHLA-TCR Analyzer - Usage Examples")
    print("=" * 50)
    
    # Create sample files first
    create_sample_input()
    
    # Run examples
    example_basic_usage()
    example_batch_processing()
    example_custom_analysis()
    
    print("\\n" + "=" * 50)
    print("Examples completed!")
    print("\\nNote: Most examples require actual PDB files to run.")
    print("Replace the example file paths with your actual PDB files.")


if __name__ == "__main__":
    main()