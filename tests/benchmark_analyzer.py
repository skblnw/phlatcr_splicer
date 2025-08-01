#!/usr/bin/env python3
"""
Benchmark script for pHLA-TCR Complex Structure Analyzer

This script tests the analyzer against known pHLA-TCR complexes
and evaluates its performance.
"""

import os
import sys
import time
import json
from collections import defaultdict
import sys
import os
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..'))
from phlatcr_splicer import pHLATCRAnalyzer


class AnalyzerBenchmark:
    """Benchmark the analyzer against known test cases."""
    
    def __init__(self):
        self.analyzer = pHLATCRAnalyzer(verbose=False)
        self.test_cases = self._load_test_cases()
        self.results = {}
    
    def _load_test_cases(self):
        """Load or define test cases with expected results."""
        # This would ideally load from a database of known structures
        # For now, we'll define some synthetic test cases
        return {
            'synthetic_1': {
                'chains': {
                    'P': {'length': 9, 'type': 'peptide', 'sequence': 'YLQPRTWYF'},
                    'B': {'length': 99, 'type': 'b2m', 'sequence': 'IQRTPK' + 'A' * 90 + 'VNHVTL'},
                    'H': {'length': 276, 'type': 'mhc_heavy', 'sequence': 'GSHSMRY' + 'G' * 260 + 'FKAFLKQ'},
                    'A': {'length': 215, 'type': 'tcr_alpha', 'sequence': 'FGXGT' + 'S' * 200 + 'WYQQKP'},
                    'T': {'length': 245, 'type': 'tcr_beta', 'sequence': 'MGIGV' + 'T' * 230 + 'FGGGT'},
                }
            },
            'synthetic_2': {
                'chains': {
                    'X': {'length': 12, 'type': 'peptide', 'sequence': 'AMETQIVAYLAQ'},
                    'Y': {'length': 100, 'type': 'b2m', 'sequence': 'SRNLTKDR' + 'F' * 85 + 'IQRTPK'},
                    'Z': {'length': 290, 'type': 'mhc_heavy', 'sequence': 'WSDRVI' + 'L' * 275 + 'TGAASC'},
                    'U': {'length': 220, 'type': 'tcr_alpha', 'sequence': 'YKFK' + 'R' * 210 + 'TLTIS'},
                    'V': {'length': 250, 'type': 'tcr_beta', 'sequence': 'SVGD' + 'E' * 240 + 'FGQGT'},
                }
            }
        }
    
    def create_test_pdb(self, test_case_name: str) -> str:
        """Create a PDB file for a test case."""
        test_case = self.test_cases[test_case_name]
        
        pdb_content = f"""HEADER    IMMUNE SYSTEM/IMMUNE SYSTEM             01-JAN-20   {test_case_name.upper()[:4]}    
TITLE     TEST CASE: {test_case_name.upper()}
"""
        
        atom_num = 1
        for chain_id, chain_data in test_case['chains'].items():
            sequence = chain_data['sequence']
            
            for i, aa_char in enumerate(sequence):
                # Convert single letter to 3-letter code (simplified)
                aa_map = {
                    'A': 'ALA', 'R': 'ARG', 'N': 'ASN', 'D': 'ASP', 'C': 'CYS',
                    'E': 'GLU', 'Q': 'GLN', 'G': 'GLY', 'H': 'HIS', 'I': 'ILE',
                    'L': 'LEU', 'K': 'LYS', 'M': 'MET', 'F': 'PHE', 'P': 'PRO',
                    'S': 'SER', 'T': 'THR', 'W': 'TRP', 'Y': 'TYR', 'V': 'VAL',
                    'X': 'ALA'  # Default for unknown
                }
                aa3 = aa_map.get(aa_char, 'ALA')
                
                # Add minimal atom records with proper PDB formatting
                x = 10.0 + i * 0.1 + ord(chain_id) * 10
                y = 10.0 + i * 0.05
                z = 10.0
                
                pdb_content += f"ATOM  {atom_num:5d}  N   {aa3} {chain_id}{i+1:4d}      {x:8.3f}{y:8.3f}{z:8.3f}  1.00 20.00           N  \n"
                atom_num += 1
                pdb_content += f"ATOM  {atom_num:5d}  CA  {aa3} {chain_id}{i+1:4d}      {x+0.5:8.3f}{y+0.5:8.3f}{z+0.5:8.3f}  1.00 20.00           C  \n"
                atom_num += 1
                pdb_content += f"ATOM  {atom_num:5d}  C   {aa3} {chain_id}{i+1:4d}      {x+1.0:8.3f}{y+1.0:8.3f}{z+1.0:8.3f}  1.00 20.00           C  \n"
                atom_num += 1
                pdb_content += f"ATOM  {atom_num:5d}  O   {aa3} {chain_id}{i+1:4d}      {x+1.5:8.3f}{y+1.5:8.3f}{z+1.5:8.3f}  1.00 20.00           O  \n"
                atom_num += 1
        
        pdb_content += "END\n"
        
        # Write to file
        filename = f"test_{test_case_name}.pdb"
        with open(filename, 'w') as f:
            f.write(pdb_content)
        
        return filename
    
    def run_benchmark(self):
        """Run benchmark on all test cases."""
        print("Running pHLA-TCR Analyzer Benchmark")
        print("=" * 50)
        
        total_start = time.time()
        
        for test_name, test_case in self.test_cases.items():
            print(f"\nTest Case: {test_name}")
            print("-" * 20)
            
            # Create test PDB
            pdb_file = self.create_test_pdb(test_name)
            
            try:
                # Time the analysis
                start_time = time.time()
                predicted = self.analyzer.analyze_pdb(pdb_file)
                analysis_time = time.time() - start_time
                
                # Evaluate results
                evaluation = self._evaluate_prediction(test_case['chains'], predicted)
                
                # Store results
                self.results[test_name] = {
                    'predicted': predicted,
                    'expected': {cid: cdata['type'] for cid, cdata in test_case['chains'].items()},
                    'evaluation': evaluation,
                    'time': analysis_time
                }
                
                # Print results
                print(f"Analysis time: {analysis_time:.3f}s")
                print(f"Accuracy: {evaluation['accuracy']:.2%}")
                print(f"Correct assignments: {evaluation['correct']}/{evaluation['total']}")
                
                for chain_id in test_case['chains']:
                    expected = test_case['chains'][chain_id]['type']
                    actual = predicted.get(chain_id, 'unknown')
                    status = "✓" if expected == actual else "✗"
                    print(f"  Chain {chain_id}: {expected} -> {actual} {status}")
                
            except Exception as e:
                print(f"Error: {e}")
                self.results[test_name] = {'error': str(e)}
            
            finally:
                # Clean up
                if os.path.exists(pdb_file):
                    os.remove(pdb_file)
        
        total_time = time.time() - total_start
        
        # Print summary
        print(f"\n{'='*50}")
        print("BENCHMARK SUMMARY")
        print(f"{'='*50}")
        print(f"Total time: {total_time:.3f}s")
        
        # Calculate overall statistics
        self._print_summary()
    
    def _evaluate_prediction(self, expected_chains, predicted):
        """Evaluate prediction against expected results."""
        correct = 0
        total = len(expected_chains)
        
        for chain_id, chain_data in expected_chains.items():
            expected_type = chain_data['type']
            predicted_type = predicted.get(chain_id, 'unknown')
            
            if expected_type == predicted_type:
                correct += 1
        
        return {
            'correct': correct,
            'total': total,
            'accuracy': correct / total if total > 0 else 0
        }
    
    def _print_summary(self):
        """Print overall benchmark summary."""
        if not self.results:
            print("No results to summarize")
            return
        
        # Calculate overall statistics
        total_correct = 0
        total_predictions = 0
        total_time = 0
        
        for test_name, result in self.results.items():
            if 'evaluation' in result:
                total_correct += result['evaluation']['correct']
                total_predictions += result['evaluation']['total']
                total_time += result['time']
        
        if total_predictions > 0:
            overall_accuracy = total_correct / total_predictions
            avg_time = total_time / len(self.results)
            
            print(f"Overall accuracy: {overall_accuracy:.2%}")
            print(f"Average analysis time: {avg_time:.3f}s")
            print(f"Total chains analyzed: {total_predictions}")
            print(f"Correct predictions: {total_correct}")
        
        # Chain type breakdown
        print("\nPer-chain-type accuracy:")
        type_stats = defaultdict(lambda: {'correct': 0, 'total': 0})
        
        for test_name, result in self.results.items():
            if 'expected' in result and 'predicted' in result:
                for chain_id, expected_type in result['expected'].items():
                    predicted_type = result['predicted'].get(chain_id, 'unknown')
                    type_stats[expected_type]['total'] += 1
                    if expected_type == predicted_type:
                        type_stats[expected_type]['correct'] += 1
        
        for chain_type, stats in type_stats.items():
            accuracy = stats['correct'] / stats['total'] if stats['total'] > 0 else 0
            print(f"  {chain_type}: {accuracy:.2%} ({stats['correct']}/{stats['total']})")
    
    def save_results(self, filename: str = "benchmark_results.json"):
        """Save benchmark results to file."""
        with open(filename, 'w') as f:
            json.dump(self.results, f, indent=2)
        print(f"\nResults saved to: {filename}")


def main():
    """Run the benchmark."""
    benchmark = AnalyzerBenchmark()
    benchmark.run_benchmark()
    benchmark.save_results()


if __name__ == "__main__":
    main()