"""
MHC-II Complex Analyzer

A specialized tool for identifying and classifying chains in MHC-II-TCR complexes.
Unlike MHC-I complexes, MHC-II complexes contain:
- MHC-II α chain (~180-200 residues)  
- MHC-II β chain (~190-210 residues)
- Longer peptide (12-25 residues)
- TCR α chain (180-230 residues)
- TCR β chain (230-290 residues)
- No β2-microglobulin

Author: AI Assistant
"""

from typing import Dict, List, Tuple, Optional, Set
from collections import defaultdict
import numpy as np
import warnings

try:
    from Bio.PDB import PDBParser, Structure, Structure, Model, Chain, Residue
    from Bio.PDB.PDBIO import PDBIO
    from Bio.SeqUtils import ProtParam
    from Bio.Seq import Seq
except ImportError:
    print("Error: BioPython is required. Install with: pip install biopython")
    raise

import os
import tempfile


class pMHCIITCRAnalyzer:
    """
    Analyzer for MHC-II-TCR complex structures.
    
    Identifies and classifies the following chain types:
    - mhc_ii_alpha: MHC-II α chain
    - mhc_ii_beta: MHC-II β chain  
    - peptide: Presented peptide antigen (longer than MHC-I)
    - tcr_alpha: T-cell receptor α chain
    - tcr_beta: T-cell receptor β chain
    """
    
    def __init__(self, verbose: bool = False):
        """
        Initialize the MHC-II analyzer.
        
        Args:
            verbose: If True, print detailed analysis information
        """
        self.verbose = verbose
        self.parser = PDBParser(QUIET=True)
        
        # Load patterns for each chain type
        self.mhc_ii_patterns = self._load_mhc_ii_patterns()
        self.tcr_patterns = self._load_tcr_patterns()
        self.peptide_patterns = self._load_peptide_patterns()
        
    def _load_mhc_ii_patterns(self) -> Dict:
        """Load MHC-II alpha and beta chain patterns and characteristics."""
        return {
            'alpha': {
                'length_range': (170, 210),  # MHC-II α chains
                'optimal_length': (180, 200),
                'highly_specific_motifs': [
                    'MGSGWV',    # Common MHC-II α N-terminal
                    'MILNKA',    # MHC-II α leader sequence variant
                    'RFDSDK',    # α1 domain pattern
                    'EPRAPW',    # α1 domain signature
                    'QFRVD',     # α1-α2 linker region
                ],
                'conserved_motifs': [
                    'HIVQR',     # α1 domain
                    'WVLVL',     # α2 domain
                    'HSLTE',     # α2 domain
                    'VTQQP',     # α2 domain
                    'CRHNY',     # α3 domain
                ],
                'n_terminal_signatures': [
                    'MGSGWV',    # Very specific α chain start
                    'MILNKA',    # Alternative start
                    'GSGWV',     # Variant without Met
                ],
                'c_terminal_patterns': [
                    'CRHNY',     # α3 domain
                    'HVTLS',     # C-terminal region
                    'GACIF',     # Transmembrane region
                ],
                'domain_boundaries': {
                    'alpha1': (1, 85),
                    'alpha2': (86, 175),
                    'alpha3': (176, 200),
                },
            },
            'beta': {
                'length_range': (180, 220),  # MHC-II β chains  
                'optimal_length': (190, 210),
                'highly_specific_motifs': [
                    'MSWKKA',    # Common MHC-II β N-terminal
                    'MIHSPS',    # MHC-II β leader sequence variant
                    'PVTVY',     # β1 domain signature
                    'NYGVV',     # β1 domain pattern
                    'WNGLR',     # β1-β2 linker
                ],
                'conserved_motifs': [
                    'FLNGQ',     # β1 domain
                    'YVGC',      # β1 domain
                    'NGTER',     # β2 domain
                    'NLLR',      # β2 domain
                    'VTKR',      # β2 domain
                ],
                'n_terminal_signatures': [
                    'MSWKKA',    # Very specific β chain start
                    'MIHSPS',    # Alternative start
                    'SWKKA',     # Variant without Met
                ],
                'c_terminal_patterns': [
                    'VTKR',      # β2 domain
                    'RPSAE',     # C-terminal region
                    'TLLV',      # Transmembrane region
                ],
                'domain_boundaries': {
                    'beta1': (1, 90),
                    'beta2': (91, 180),
                    'beta3': (181, 210),
                },
            }
        }
    
    def _load_tcr_patterns(self) -> Dict:
        """Load TCR alpha and beta chain patterns (similar to MHC-I analyzer)."""
        return {
            'alpha': {
                'length_range': (180, 230),  # Alpha chains are typically shorter
                'optimal_length': (190, 220),
                'highly_specific_motifs': [
                    'WYQQFP',    # Very specific TCR alpha pattern
                    'WYQQ',      # Alpha-specific framework
                    'TNDYIT',    # Alpha V-region pattern
                    'QPISMDS',   # Alpha N-terminal pattern
                    'LAKTTQ',    # Specific alpha start
                ],
                'conserved_motifs': [
                    'FGXGT',    # CDR3 C-terminal
                    'WYQQKP',   # Framework regions
                    'YKFK',     # Framework 1
                    'TLTIS',    # Framework 3
                    'QLLE',     # Common N-terminal for alpha
                    'SPQFL',    # Alpha-specific pattern
                    'DSYEG',    # Alpha framework
                ],
                'n_terminal_patterns': ['QLLE', 'QVQL', 'EVQL', 'LAKTTQ'],  # Alpha N-terminal
                'c_terminal_patterns': ['FGXGT', 'FGAGT', 'FGQGT', 'FPSSP'],
            },
            'beta': {
                'length_range': (230, 290),  # Beta chains are typically longer
                'optimal_length': (240, 280),
                'highly_specific_motifs': [
                    'VKVTQS',    # Very specific beta N-terminal
                    'SSRYLV',    # Beta framework pattern
                    'LIYFSYD',   # Beta CDR2 region
                    'RTGEKV',    # Beta V-region
                    'YRQDPG',    # Beta framework
                ],
                'conserved_motifs': [
                    'FGGGT',    # CDR3 C-terminal
                    'WYQQKP',   # Framework regions
                    'MGIGV',    # Framework 1
                    'SVGD',     # Framework 2
                    'GITQ',     # Common N-terminal for beta
                    'SPKYL',    # Beta-specific pattern
                ],
                'n_terminal_patterns': ['GITQ', 'MGIT', 'GVTQ', 'VKVTQ'],  # Beta N-terminal
                'c_terminal_patterns': ['FGGGT', 'FGQGT', 'FGPGT', 'WTQDRA'],
            }
        }
    
    def _load_peptide_patterns(self) -> Dict:
        """Load peptide patterns for MHC-II (longer peptides than MHC-I)."""
        return {
            'length_range': (12, 25),  # MHC-II peptides are longer
            'optimal_length': (14, 20),  # Most common range
            'common_residues': ['A', 'L', 'V', 'I', 'F', 'Y', 'W'],  # Hydrophobic anchor residues
            'molecular_weight': (1200, 3000),  # Da
        }
    
    def analyze_pdb(self, pdb_file: str) -> Dict[str, str]:
        """
        Analyze a PDB file and identify MHC-II-TCR complex chain types.
        
        Args:
            pdb_file: Path to the PDB file
            
        Returns:
            Dictionary mapping chain IDs to chain types
        """
        if not os.path.exists(pdb_file):
            raise FileNotFoundError(f"PDB file not found: {pdb_file}")
        
        if self.verbose:
            print(f"Analyzing MHC-II complex PDB file: {pdb_file}")
        
        # Parse the structure
        structure = self.parser.get_structure('complex', pdb_file)
        
        # Extract chains and their sequences
        chain_info = self._extract_chain_info(structure)
        
        if self.verbose:
            print(f"Found {len(chain_info)} chains")
            for chain_id, info in chain_info.items():
                print(f"  Chain {chain_id}: {len(info['sequence'])} residues")
        
        # Scenario detection based on header presence
        header_present = self._has_header_section(pdb_file)
        header_preassignments = self._extract_header_chain_annotations(pdb_file) if header_present else {}
        
        # If header exists, prefer REMARK-based grouping; else assume single complex
        if header_present:
            # Count REMARK 350 chain-group lines; if present, trust REMARK grouping even if single set
            remark_group_count = self._count_remark350_groups(pdb_file)
            # Try REMARK-based grouping first (REMARK 350 chain groups)
            complexes = self._detect_complexes_from_remarks(pdb_file, chain_info)
            # Only fall back to heuristics if there were no REMARK groups at all
            if remark_group_count == 0:
                complexes = self._detect_complexes(chain_info, structure)
        else:
            complexes = [chain_info]
        
        if self.verbose and len(complexes) > 1:
            print(f"\nDetected {len(complexes)} complexes:")
            for i, complex_chains in enumerate(complexes, 1):
                chain_ids = sorted(complex_chains.keys())
                print(f"  Complex {i}: Chains {', '.join(chain_ids)}")
        
        # Classify chains within each complex
        if len(complexes) > 1:
            chain_assignments = self._classify_multiple_complexes(complexes)
            # Apply header-based preassignments if available (with sanity checks)
            if header_preassignments:
                chain_assignments = self._apply_header_preassignments(chain_assignments, complexes, chain_info, header_preassignments)
        else:
            # Single complex - use original logic
            chain_assignments = self._classify_chains(chain_info)
            # Apply header-based preassignments if available
            if header_preassignments:
                single_complexes = [chain_info]
                chain_assignments = self._apply_header_preassignments(chain_assignments, single_complexes, chain_info, header_preassignments)
        
        # If no header present (likely predicted single complex), validate and warn on failure
        if not header_present:
            self._validate_single_complex_expectations(chain_assignments)
        
        if self.verbose:
            print("\nFinal chain assignments:")
            # Sort by chain ID for consistent output
            for chain_id in sorted(chain_assignments.keys()):
                print(f"  Chain {chain_id}: {chain_assignments[chain_id]}")
        
        # Emit warnings if required components are missing (per complex)
        self._emit_component_warnings(chain_assignments, header_preassignments)
        
        return chain_assignments

    def _apply_header_preassignments(self, chain_assignments: Dict[str, str], complexes: List[Dict], chain_info: Dict, header_preassignments: Dict[str, str]) -> Dict[str, str]:
        """Apply header annotations cautiously: only when length-consistent and helpful.
        - Do not override a non-unknown assignment unless types match and are consistent by length
        - For peptides, ensure chain length <= 30
        - For MHC-II alpha/beta, ensure lengths are in broad expected ranges
        - For TCR alpha/beta, ensure lengths are roughly in expected ranges
        """
        # Build reverse map chain -> complex index
        chain_to_complex = {}
        for idx, comp in enumerate(complexes, start=1):
            for cid in comp.keys():
                chain_to_complex[cid] = idx
        
        def length_ok(ctype: str, length: int) -> bool:
            if ctype == 'peptide':
                return length <= 30
            if ctype == 'mhc_ii_alpha':
                return 160 <= length <= 230
            if ctype == 'mhc_ii_beta':
                return 160 <= length <= 240
            if ctype == 'tcr_alpha':
                # Accept full (180-240) and V-domain-only (100-160)
                return (170 <= length <= 240) or (100 <= length <= 160)
            if ctype == 'tcr_beta':
                # Accept full (230-320) and V-domain-only (100-170)
                return (230 <= length <= 320) or (100 <= length <= 170)
            return True
        
        updated = dict(chain_assignments)
        for cid, annotated_type in header_preassignments.items():
            if cid not in updated or cid not in chain_info:
                continue
            current = updated[cid]
            length = chain_info[cid]['length']
            if not length_ok(annotated_type, length):
                continue
            base_current = current.split('_complex')[0]
            # For MHC-II alpha/beta and peptide, trust header over scoring
            if annotated_type in ('mhc_ii_alpha', 'mhc_ii_beta', 'peptide'):
                comp_idx = chain_to_complex.get(cid)
                if comp_idx is not None and len(complexes) > 1:
                    updated[cid] = f"{annotated_type}_complex{comp_idx}"
                else:
                    updated[cid] = annotated_type
                continue
            # For TCR types, apply only if current is unknown or matches annotated type
            if base_current == 'unknown' or base_current == annotated_type:
                comp_idx = chain_to_complex.get(cid)
                if comp_idx is not None and len(complexes) > 1:
                    updated[cid] = f"{annotated_type}_complex{comp_idx}"
                else:
                    updated[cid] = annotated_type
        return updated

    def _detect_complexes_from_remarks(self, pdb_file: str, chain_info: Dict) -> List[Dict]:
        """
        Use REMARK 350 lines to pre-group chains into complexes when available.
        
        In multi-complex entries, lines like "APPLY THE FOLLOWING TO CHAINS: A, B, ..."
        hint at biological assembly groupings that often correspond to a single complex.
        """
        try:
            groups: List[List[str]] = []
            with open(pdb_file, 'r') as fh:
                for line in fh:
                    if line.startswith('REMARK 350') and 'APPLY THE FOLLOWING TO CHAINS:' in line:
                        parts = line.split('APPLY THE FOLLOWING TO CHAINS:')
                        if len(parts) < 2:
                            continue
                        tail = parts[1].strip()
                        tail = tail.split('  ')[0].strip()
                        raw_tokens = [t.strip() for t in tail.split(',') if t.strip()]
                        chains = [tok[0] for tok in raw_tokens if len(tok) > 0]
                        if len(chains) >= 4:
                            groups.append(chains)
            # Deduplicate identical groups
            unique_groups: List[List[str]] = []
            seen = set()
            for grp in groups:
                key = tuple(sorted(grp))
                if key not in seen:
                    seen.add(key)
                    unique_groups.append(grp)
            # Remove subset groups (keep only maximal sets)
            maximal_groups: List[List[str]] = []
            unique_sets = [set(g) for g in unique_groups]
            for i, s in enumerate(unique_sets):
                is_subset = False
                for j, t in enumerate(unique_sets):
                    if i != j and s.issubset(t):
                        is_subset = True
                        break
                if not is_subset:
                    maximal_groups.append(sorted(list(s)))
            # If overlapping groups remain (share chains), treat as single complex
            for i in range(len(maximal_groups)):
                for j in range(i+1, len(maximal_groups)):
                    if set(maximal_groups[i]).intersection(set(maximal_groups[j])):
                        return [chain_info]
            # Map groups to chain_info
            if len(maximal_groups) >= 2:
                complexes: List[Dict] = []
                for grp in maximal_groups:
                    complex_chains: Dict = {}
                    for cid in grp:
                        if cid in chain_info:
                            complex_chains[cid] = chain_info[cid]
                    if complex_chains:
                        complexes.append(complex_chains)
                if len([c for c in complexes if len(c) >= 2]) >= 2:
                    return complexes
            return [chain_info]
        except Exception:
            return [chain_info]

    def _count_remark350_groups(self, pdb_file: str) -> int:
        """Count the number of REMARK 350 lines that apply to chains."""
        count = 0
        try:
            with open(pdb_file, 'r') as fh:
                for line in fh:
                    if line.startswith('REMARK 350') and 'APPLY THE FOLLOWING TO CHAINS:' in line:
                        count += 1
        except Exception:
            return 0
        return count

    def _has_header_section(self, pdb_file: str) -> bool:
        """Detect whether the file has a PDB header section before ATOM/HETATM lines."""
        try:
            with open(pdb_file, 'r') as fh:
                for line in fh:
                    rec = line[:6]
                    if rec.startswith('ATOM') or rec.startswith('HETATM'):
                        return False
                    if rec.startswith('HEADER') or rec.startswith('TITLE') or rec.startswith('COMPND') or rec.startswith('REMARK') or rec.startswith('DBREF') or rec.startswith('SEQRES'):
                        return True
        except Exception:
            return False
        return False

    def _extract_header_chain_annotations(self, pdb_file: str) -> Dict[str, str]:
        """
        Extract chain type hints from PDB header (COMPND, DBREF, etc.).
        Returns a mapping of chain_id -> chain_type.
        """
        annotations: Dict[str, str] = {}
        try:
            current_mol = None
            mol_to_chains: Dict[str, List[str]] = {}
            mol_to_desc: Dict[str, str] = {}
            with open(pdb_file, 'r') as fh:
                for line in fh:
                    if line.startswith('ATOM') or line.startswith('HETATM'):
                        break
                    if line.startswith('COMPND'):
                        content = line[10:].strip()
                        # Parse MOL_ID
                        if 'MOL_ID:' in content:
                            parts = content.split('MOL_ID:')
                            rest = parts[1].strip()
                            current_mol = rest.split(';')[0].strip()
                        if 'MOLECULE:' in content and current_mol is not None:
                            desc = content.split('MOLECULE:')[1].strip()
                            desc = desc.split(';')[0].strip()
                            mol_to_desc[current_mol] = mol_to_desc.get(current_mol, '') + ' ' + desc
                        if 'CHAIN:' in content and current_mol is not None:
                            chains_field = content.split('CHAIN:')[1].strip()
                            chains_field = chains_field.split(';')[0]
                            chains = [c.strip() for c in chains_field.replace('.', '').split(',') if c.strip()]
                            mol_to_chains.setdefault(current_mol, [])
                            for c in chains:
                                if c and c[0].isalpha():
                                    mol_to_chains[current_mol].append(c[0])
                    elif line.startswith('DBREF'):
                        # DBREF lines include chain and DB identifiers and sequence range
                        try:
                            tokens = line.split()
                            # tokens: DBREF, pdb, chain, start, end, db, acc, dbname, dbstart, dbend
                            chain_id = tokens[2]
                            start = int(tokens[3]) if tokens[3].lstrip('-').isdigit() else None
                            end = int(tokens[4]) if tokens[4].lstrip('-').isdigit() else None
                            # Map by DB name for general types
                            db_id = ' '.join(tokens[5:]).upper()
                            mapped = self._map_dbref_to_type(db_id)
                            if mapped:
                                annotations[chain_id] = mapped
                            # If range is short (<=30), treat as peptide regardless of db_id
                            if start is not None and end is not None:
                                seg_len = abs(end - start) + 1
                                if seg_len <= 30:
                                    annotations[chain_id] = 'peptide'
                        except Exception:
                            # Fallback to previous simple parsing
                            chain_id = line[12].strip()
                            db_id = line[42:].strip().upper()
                            mapped = self._map_dbref_to_type(db_id)
                            if mapped:
                                annotations[chain_id] = mapped
            # Map COMPND molecules to types
            for mol_id, chains in mol_to_chains.items():
                desc = mol_to_desc.get(mol_id, '').upper()
                mapped = self._map_molecule_desc_to_type(desc)
                if mapped:
                    for cid in chains:
                        annotations.setdefault(cid, mapped)
                # Heuristic: if molecule description suggests antigen/protein (e.g., TENASCIN) and not HLA/TCR,
                # mark as peptide; will be validated by length later when applying
                if ('TENASCIN' in desc or 'FRAGMENT' in desc or 'RESIDUES' in desc) and not any(k in desc for k in ['HLA', 'MHC', 'TCR']):
                    for cid in chains:
                        annotations.setdefault(cid, 'peptide')
        except Exception:
            return annotations
        return annotations

    def _map_molecule_desc_to_type(self, desc: str) -> Optional[str]:
        """Map COMPND MOLECULE description to chain type for MHC-II."""
        if not desc:
            return None
        # Accept generic TCR descriptions as well as explicit T-CELL RECEPTOR
        if ('T-CELL' in desc or 'T CELL' in desc or 'TCR' in desc) and ('RECEPTOR' in desc or 'TCR' in desc or 'CHAIN' in desc):
            if 'ALPHA' in desc:
                return 'tcr_alpha'
            if 'BETA' in desc:
                return 'tcr_beta'
            # Heuristic: "TCR A CHAIN" / "TCR B CHAIN"
            if ' A CHAIN' in desc or desc.endswith(' TCR A') or ' TCR A ' in desc:
                return 'tcr_alpha'
            if ' B CHAIN' in desc or desc.endswith(' TCR B') or ' TCR B ' in desc:
                return 'tcr_beta'
        if 'HLA-DR' in desc or 'HLA-DQ' in desc or 'HLA-DP' in desc or 'MHC CLASS II' in desc or 'CLASS II' in desc:
            # Need to distinguish α vs β later; COMPND may list separately
            if 'ALPHA' in desc:
                return 'mhc_ii_alpha'
            if 'BETA' in desc:
                return 'mhc_ii_beta'
        if 'PEPTIDE' in desc or 'ANTIGEN' in desc or 'EPITOPE' in desc:
            return 'peptide'
        return None

    def _map_dbref_to_type(self, db_id: str) -> Optional[str]:
        """Map DBREF db_id (e.g., UNP name) to chain type for MHC-II."""
        if not db_id:
            return None
        if 'TRAC' in db_id or 'TRA' in db_id or 'TRAV' in db_id or 'TCRA' in db_id:
            return 'tcr_alpha'
        if 'TRBC' in db_id or 'TRB' in db_id or 'TRBV' in db_id or 'TCRB' in db_id:
            return 'tcr_beta'
        if 'HLA-DRA' in db_id or 'DRA_' in db_id or 'HLA-DQA' in db_id or 'HLA-DPA' in db_id or 'HLA-D' in db_id and 'A' in db_id:
            return 'mhc_ii_alpha'
        if 'HLA-DRB' in db_id or 'DRB_' in db_id or 'HLA-DQB' in db_id or 'HLA-DPB' in db_id or 'HLA-D' in db_id and 'B' in db_id:
            return 'mhc_ii_beta'
        return None

    def _validate_single_complex_expectations(self, assignments: Dict[str, str]) -> None:
        """
        Validate minimal expectations for a single MHC-II complex when no header is present.
        Raise a RuntimeError with a warning message if expectations are not met.
        """
        # Minimal expected: at least one mhc_ii_alpha, one mhc_ii_beta, one peptide, one tcr_alpha, one tcr_beta
        required = ['mhc_ii_alpha', 'mhc_ii_beta', 'peptide', 'tcr_alpha', 'tcr_beta']
        present = {t: False for t in required}
        for v in assignments.values():
            base = v.split('_complex')[0]
            if base in present:
                present[base] = True
        if not all(present.values()):
            missing = [k for k, ok in present.items() if not ok]
            raise RuntimeError(
                f"No header found; expected a single pMHC-II-TCR complex, but missing required components: {', '.join(missing)}"
            )
    
    def _emit_component_warnings(self, assignments: Dict[str, str], header_preassignments: Optional[Dict[str, str]] = None) -> None:
        """
        Issue warnings if any complex is missing required components.
        This runs regardless of header presence, but does not raise.
        """
        # Determine if TCR is expected based on header annotations (if provided)
        expect_tcr = False
        if header_preassignments:
            expect_tcr = any(v in ('tcr_alpha', 'tcr_beta') or v.startswith('tcr_') for v in header_preassignments.values())
        # Always require MHC-II α/β and peptide; require TCR only when header indicates presence
        required = ['mhc_ii_alpha', 'mhc_ii_beta', 'peptide'] + (['tcr_alpha', 'tcr_beta'] if expect_tcr else [])
        # Group by complex index parsed from suffix _complexN (default to 1)
        complexes: Dict[int, Dict[str, bool]] = defaultdict(lambda: {t: False for t in required})
        for chain_id, ctype in assignments.items():
            base = ctype.split('_complex')[0]
            comp_idx = 1
            if '_complex' in ctype:
                try:
                    comp_idx = int(ctype.split('_complex')[-1])
                except Exception:
                    comp_idx = 1
            if base in complexes[comp_idx]:
                complexes[comp_idx][base] = True
        # Emit warnings for any missing components
        missing_summary = {}
        for idx, present in complexes.items():
            missing = [t for t, ok in present.items() if not ok]
            if missing:
                missing_summary[idx] = missing
                warnings.warn(
                    f"MHC-II complex {idx} missing required component(s): {', '.join(missing)}",
                    UserWarning
                )
        # If anything missing anywhere, stop and raise an error
        if missing_summary:
            parts = [f"complex {idx}: {', '.join(miss)}" for idx, miss in sorted(missing_summary.items())]
            raise RuntimeError(
                "Missing required component(s) in MHC-II analysis (" + "; ".join(parts) + ")"
            )
    
    def _extract_chain_info(self, structure: Structure) -> Dict:
        """
        Extract information about each chain in the structure.
        
        Args:
            structure: BioPython structure object
            
        Returns:
            Dictionary with chain information
        """
        chain_info = {}
        
        for model in structure:
            for chain in model:
                chain_id = chain.id
                sequence = self._extract_sequence(chain)
                
                if len(sequence) > 0:  # Only include chains with sequence
                    # Calculate molecular properties
                    properties = self._calculate_properties(sequence)
                    
                    chain_info[chain_id] = {
                        'sequence': sequence,
                        'length': len(sequence),
                        'properties': properties
                    }
        
        return chain_info
    
    def _extract_sequence(self, chain) -> str:
        """Extract amino acid sequence from a chain."""
        sequence = []
        
        for residue in chain:
            if residue.has_id('CA'):  # Only standard amino acids
                res_name = residue.get_resname().strip()
                aa = self._three_to_one(res_name)
                if aa != 'X':  # Valid amino acid
                    sequence.append(aa)
        
        return ''.join(sequence)
    
    def _three_to_one(self, three_letter: str) -> str:
        """Convert three-letter amino acid code to one-letter."""
        aa_map = {
            'ALA': 'A', 'CYS': 'C', 'ASP': 'D', 'GLU': 'E', 'PHE': 'F',
            'GLY': 'G', 'HIS': 'H', 'ILE': 'I', 'LYS': 'K', 'LEU': 'L',
            'MET': 'M', 'ASN': 'N', 'PRO': 'P', 'GLN': 'Q', 'ARG': 'R',
            'SER': 'S', 'THR': 'T', 'VAL': 'V', 'TRP': 'W', 'TYR': 'Y'
        }
        return aa_map.get(three_letter, 'X')
    
    def _calculate_properties(self, sequence: str) -> Dict:
        """Calculate molecular properties of a sequence."""
        try:
            protein = ProtParam.ProteinAnalysis(sequence)
            return {
                'molecular_weight': protein.molecular_weight(),
                'gravy': protein.gravy(),  # Hydrophobicity
                'instability_index': protein.instability_index(),
                'charge': sum(1 for aa in sequence if aa in 'KR') - sum(1 for aa in sequence if aa in 'DE')
            }
        except Exception:
            return {
                'molecular_weight': len(sequence) * 110,  # Average AA weight
                'gravy': 0.0,
                'instability_index': 40.0,
                'charge': 0
            }
    
    def _detect_complexes(self, chain_info: Dict, structure) -> List[Dict]:
        """
        Detect multiple MHC-II-TCR complexes in the structure.
        Uses similar logic to MHC-I analyzer but adapted for MHC-II.
        """
        # First, try chain counting approach
        complexes_by_counting = self._detect_complexes_by_counting(chain_info)
        
        if self.verbose and len(complexes_by_counting) > 1:
            print(f"Chain counting suggests {len(complexes_by_counting)} complexes")
        
        # If counting detects multiple complexes, use that result
        if len(complexes_by_counting) > 1:
            return complexes_by_counting
        
        # If counting suggests single complex but we have many chains, try spatial clustering
        if len(chain_info) >= 8:  # Only try spatial if we have many chains
            complexes_by_spatial = self._detect_complexes_by_spatial_clustering(chain_info, structure)
            
            if self.verbose and len(complexes_by_spatial) > 1:
                print(f"Spatial clustering suggests {len(complexes_by_spatial)} complexes")
            
            return complexes_by_spatial
        
        # Default to single complex
        return complexes_by_counting
    
    def _detect_complexes_by_counting(self, chain_info: Dict) -> List[Dict]:
        """
        Detect complexes by counting potential chain types for MHC-II.
        Enhanced algorithm that understands MHC-II α/β pairing requirements.
        """
        # Score all chains for all types
        potential_types = {
            'peptide': [],
            'mhc_ii_alpha': [],
            'mhc_ii_beta': [],
            'tcr_alpha': [],
            'tcr_beta': []
        }
        
        all_scores = {}
        
        for chain_id, info in chain_info.items():
            seq = info['sequence']
            length = len(seq)
            props = info['properties']
            
            # Score for each type
            scores = {
                'peptide': self._score_peptide(seq, length, props),
                'mhc_ii_alpha': self._score_mhc_ii_alpha(seq, length),
                'mhc_ii_beta': self._score_mhc_ii_beta(seq, length),
                'tcr_alpha': self._score_tcr_alpha(seq, length),
                'tcr_beta': self._score_tcr_beta(seq, length)
            }
            
            all_scores[chain_id] = scores
            
            # Add to potential types if score is reasonable
            for chain_type, score in scores.items():
                if score > 0.3:  # Threshold for being a potential candidate
                    potential_types[chain_type].append((chain_id, score))
        
        # Sort by score for each type
        for chain_type in potential_types:
            potential_types[chain_type].sort(key=lambda x: x[1], reverse=True)
        
        # Enhanced MHC-II complex detection
        # Each MHC-II complex ideally has: 1 α + 1 β + 1 peptide + 1 TCR α + 1 TCR β
        
        # Count clear peptides (most reliable for complex count)
        clear_peptides = [cid for cid, score in potential_types['peptide'] if score > 0.8]
        
        # Estimate number of complexes from peptides; fallback to 1
        estimated_complexes = len(clear_peptides) if clear_peptides else 1
        
        if self.verbose:
            print(f"MHC-II complex estimation: {estimated_complexes} complexes based on {len(clear_peptides)} clear peptides")
        
        if estimated_complexes == 1:
            return [chain_info]
        
        # Build complexes using sophisticated pairing
        complexes = self._build_mhc_ii_complexes(chain_info, all_scores, estimated_complexes)
        
        return complexes if len(complexes) > 1 else [chain_info]
    
    def _build_mhc_ii_complexes(self, chain_info: Dict, all_scores: Dict, num_complexes: int) -> List[Dict]:
        """
        Build MHC-II complexes using sophisticated pairing algorithm.
        """
        complexes = []
        assigned_chains = set()
        
        # Get candidates for each type
        peptide_candidates = [(cid, all_scores[cid]['peptide']) for cid in chain_info.keys() 
                             if all_scores[cid]['peptide'] > 0.6]
        peptide_candidates.sort(key=lambda x: x[1], reverse=True)
        
        # For each estimated complex, try to build a complete set
        for i in range(num_complexes):
            complex_chains = {}
            
            # 1. Assign peptide (highest scoring available)
            peptide_assigned = False
            for chain_id, score in peptide_candidates:
                if chain_id not in assigned_chains:
                    complex_chains[chain_id] = chain_info[chain_id]
                    assigned_chains.add(chain_id)
                    peptide_assigned = True
                    if self.verbose:
                        print(f"  Complex {i+1}: Assigned peptide {chain_id} (score: {score:.3f})")
                    break
            
            # 2. Find best MHC-II α chain (that hasn't been assigned)
            best_alpha = None
            best_alpha_score = 0
            for chain_id in chain_info.keys():
                if (chain_id not in assigned_chains and 
                    all_scores[chain_id]['mhc_ii_alpha'] > best_alpha_score and
                    all_scores[chain_id]['mhc_ii_alpha'] > 0.5):
                    best_alpha = chain_id
                    best_alpha_score = all_scores[chain_id]['mhc_ii_alpha']
            
            if best_alpha:
                complex_chains[best_alpha] = chain_info[best_alpha]
                assigned_chains.add(best_alpha)
                if self.verbose:
                    print(f"  Complex {i+1}: Assigned MHC-II α {best_alpha} (score: {best_alpha_score:.3f})")
            
            # 3. Find best MHC-II β chain (prefer chains that score poorly for α but reasonably for β)
            best_beta = None
            best_beta_score = 0
            for chain_id in chain_info.keys():
                if chain_id not in assigned_chains:
                    alpha_score = all_scores[chain_id]['mhc_ii_alpha']
                    beta_score = all_scores[chain_id]['mhc_ii_beta']
                    tcr_alpha_score = all_scores[chain_id]['tcr_alpha']
                    tcr_beta_score = all_scores[chain_id]['tcr_beta']
                    length = len(chain_info[chain_id]['sequence'])
                    
                    # Enhanced MHC-II β selection logic
                    # Prefer chains that:
                    # 1. Are in extended MHC-II length range (160-220)
                    # 2. Score lower for α than for β, OR have reasonable β score, OR fit length profile
                    # 3. Don't score highly for TCR
                    # 4. Haven't been assigned yet
                    is_mhc_length = 160 <= length <= 220
                    prefers_beta = beta_score >= alpha_score
                    reasonable_beta = beta_score > 0.1
                    not_strong_tcr = max(tcr_alpha_score, tcr_beta_score) < 0.7
                    
                    # Calculate composite score for β chain selection
                    composite_beta_score = beta_score
                    if is_mhc_length:
                        composite_beta_score += 0.2
                    if prefers_beta:
                        composite_beta_score += 0.1
                    if not_strong_tcr:
                        composite_beta_score += 0.1
                    
                    if (is_mhc_length and (prefers_beta or reasonable_beta) and 
                        not_strong_tcr and composite_beta_score > best_beta_score):
                        best_beta = chain_id
                        best_beta_score = composite_beta_score
            
            if best_beta:
                complex_chains[best_beta] = chain_info[best_beta]
                assigned_chains.add(best_beta)
                if self.verbose:
                    print(f"  Complex {i+1}: Assigned MHC-II β {best_beta} (score: {best_beta_score:.3f})")
            
            # 4. Assign TCR chains (highest scoring available)
            # TCR α
            best_tcr_alpha = None
            best_tcr_alpha_score = 0
            for chain_id in chain_info.keys():
                if (chain_id not in assigned_chains and 
                    all_scores[chain_id]['tcr_alpha'] > best_tcr_alpha_score and
                    all_scores[chain_id]['tcr_alpha'] > 0.5):
                    best_tcr_alpha = chain_id
                    best_tcr_alpha_score = all_scores[chain_id]['tcr_alpha']
            
            if best_tcr_alpha:
                complex_chains[best_tcr_alpha] = chain_info[best_tcr_alpha]
                assigned_chains.add(best_tcr_alpha)
                if self.verbose:
                    print(f"  Complex {i+1}: Assigned TCR α {best_tcr_alpha} (score: {best_tcr_alpha_score:.3f})")
            
            # TCR β
            best_tcr_beta = None
            best_tcr_beta_score = 0
            for chain_id in chain_info.keys():
                if (chain_id not in assigned_chains and 
                    all_scores[chain_id]['tcr_beta'] > best_tcr_beta_score and
                    all_scores[chain_id]['tcr_beta'] > 0.5):
                    best_tcr_beta = chain_id
                    best_tcr_beta_score = all_scores[chain_id]['tcr_beta']
            
            if best_tcr_beta:
                complex_chains[best_tcr_beta] = chain_info[best_tcr_beta]
                assigned_chains.add(best_tcr_beta)
                if self.verbose:
                    print(f"  Complex {i+1}: Assigned TCR β {best_tcr_beta} (score: {best_tcr_beta_score:.3f})")
            
            if complex_chains:  # Only add if we found at least some chains
                complexes.append(complex_chains)
        
        # Handle any remaining unassigned chains
        remaining_chains = {cid: info for cid, info in chain_info.items() if cid not in assigned_chains}
        if remaining_chains:
            if complexes:
                # Add to existing complexes based on best fit
                for chain_id, info in remaining_chains.items():
                    best_complex_idx = 0
                    # Add to the complex with the fewest chains
                    for i, complex_chains in enumerate(complexes):
                        if len(complex_chains) < len(complexes[best_complex_idx]):
                            best_complex_idx = i
                    complexes[best_complex_idx][chain_id] = info
                    if self.verbose:
                        print(f"  Added remaining chain {chain_id} to complex {best_complex_idx + 1}")
            else:
                # No complexes were built, return all as single complex
                complexes = [chain_info]
        
        return complexes
    
    def _detect_complexes_by_spatial_clustering(self, chain_info: Dict, structure) -> List[Dict]:
        """
        Detect complexes by spatial clustering of chain centers.
        (Similar to MHC-I analyzer)
        """
        try:
            # Calculate center of mass for each chain
            chain_centers = {}
            
            for model in structure:
                for chain in model:
                    chain_id = chain.id
                    if chain_id in chain_info:
                        # Calculate center of mass
                        coords = []
                        for residue in chain:
                            if residue.has_id('CA'):  # Use alpha carbon
                                ca_atom = residue['CA']
                                coords.append(ca_atom.coord)
                        
                        if coords:
                            center = np.mean(coords, axis=0)
                            chain_centers[chain_id] = center
            
            if len(chain_centers) < 2:
                return [chain_info]
            
            # Calculate pairwise distances
            chain_ids = list(chain_centers.keys())
            distances = {}
            
            for i, chain1 in enumerate(chain_ids):
                for j, chain2 in enumerate(chain_ids[i+1:], i+1):
                    dist = np.linalg.norm(chain_centers[chain1] - chain_centers[chain2])
                    distances[(chain1, chain2)] = dist
            
            # Simple clustering: chains within 50Å are in the same complex
            cluster_threshold = 50.0
            clusters = []
            assigned = set()
            
            for chain_id in chain_ids:
                if chain_id in assigned:
                    continue
                
                cluster = {chain_id}
                to_check = [chain_id]
                
                while to_check:
                    current = to_check.pop()
                    for other_chain in chain_ids:
                        if other_chain in cluster or other_chain in assigned:
                            continue
                        
                        # Check distance
                        pair = tuple(sorted([current, other_chain]))
                        if pair in distances and distances[pair] < cluster_threshold:
                            cluster.add(other_chain)
                            to_check.append(other_chain)
                
                clusters.append(cluster)
                assigned.update(cluster)
            
            # Convert clusters to complex format
            complexes = []
            for cluster in clusters:
                complex_chains = {cid: chain_info[cid] for cid in cluster}
                complexes.append(complex_chains)
            
            return complexes if len(complexes) > 1 else [chain_info]
        
        except Exception as e:
            if self.verbose:
                print(f"Spatial clustering failed: {e}, falling back to single complex")
            return [chain_info]
    
    def _classify_chains(self, chain_info: Dict) -> Dict[str, str]:
        """
        Classify chains in a single MHC-II complex.
        """
        assignments = {}
        
        # Two-pass classification system
        # Pass 1: Calculate all scores
        all_scores = {}
        for chain_id, info in chain_info.items():
            seq = info['sequence']
            length = len(seq)
            props = info['properties']
            
            scores = {
                'peptide': self._score_peptide(seq, length, props),
                'mhc_ii_alpha': self._score_mhc_ii_alpha(seq, length),
                'mhc_ii_beta': self._score_mhc_ii_beta(seq, length),
                'tcr_alpha': self._score_tcr_alpha(seq, length),
                'tcr_beta': self._score_tcr_beta(seq, length)
            }
            all_scores[chain_id] = scores
        
        # Pass 2: Assign chains based on best scores, avoiding conflicts
        used_unique_types = set()
        
        # Sort chains by length (shortest first for peptide priority)
        sorted_chains = sorted(chain_info.items(), key=lambda x: x[1]['length'])
        
        for chain_id, info in sorted_chains:
            scores = all_scores[chain_id]
            
            # Find best scoring type that hasn't been used (for unique types)
            best_type = max(scores.items(), key=lambda x: x[1])
            best_score = best_type[1]
            best_chain_type = best_type[0]
            
            # Check if this is a unique type already assigned
            # For MHC-II complexes, each of these should appear at most once per complex
            unique_types = {'peptide', 'mhc_ii_alpha', 'mhc_ii_beta', 'tcr_alpha', 'tcr_beta'}
            
            if best_chain_type in unique_types and best_chain_type in used_unique_types:
                # Find next best type that's not unique or not used
                sorted_scores = sorted(scores.items(), key=lambda x: x[1], reverse=True)
                assigned = False
                for chain_type, score in sorted_scores:
                    # Slightly relax threshold for TCR types to avoid duplicate alpha assignment
                    min_thresh = 0.25 if chain_type in {'tcr_alpha', 'tcr_beta'} else 0.3
                    if score > min_thresh:
                        if chain_type not in unique_types or chain_type not in used_unique_types:
                            assignments[chain_id] = chain_type
                            if chain_type in unique_types:
                                used_unique_types.add(chain_type)
                            assigned = True
                            break
                
                if not assigned:
                    assignments[chain_id] = 'unknown'
            else:
                # Can assign the best type
                if best_score > 0.3:
                    assignments[chain_id] = best_chain_type
                    if best_chain_type in unique_types:
                        used_unique_types.add(best_chain_type)
                else:
                    assignments[chain_id] = 'unknown'
        
        # Post-process to ensure one TCR alpha and one TCR beta per complex when possible
        try:
            tcr_alpha_present = any(ct == 'tcr_alpha' for ct in assignments.values())
            tcr_beta_present = any(ct == 'tcr_beta' for ct in assignments.values())
            
            # If beta missing, try to assign one unknown chain as beta based on length/score
            if not tcr_beta_present:
                beta_candidates = []
                for cid, info in chain_info.items():
                    if assignments.get(cid, 'unknown') == 'unknown':
                        length = info.get('length', 0)
                        if 215 <= length <= 300:
                            seq = info.get('sequence', '')
                            beta_candidates.append((cid, self._score_tcr_beta(seq, length)))
                if beta_candidates:
                    # Pick highest scoring candidate
                    beta_candidates.sort(key=lambda x: x[1], reverse=True)
                    best_cid, best_score = beta_candidates[0]
                    if best_score >= 0.2:
                        assignments[best_cid] = 'tcr_beta'
                    else:
                        # Length-based fallback: choose the longest unknown as beta
                        longest_unknown = max(
                            ((cid, info.get('length', 0)) for cid, info in chain_info.items() if assignments.get(cid, 'unknown') == 'unknown'),
                            key=lambda x: x[1],
                            default=None
                        )
                        if longest_unknown and longest_unknown[1] >= 220:
                            assignments[longest_unknown[0]] = 'tcr_beta'
            
            # If alpha missing, try to assign one unknown chain as alpha based on length/score
            if not tcr_alpha_present:
                alpha_candidates = []
                for cid, info in chain_info.items():
                    if assignments.get(cid, 'unknown') == 'unknown':
                        length = info.get('length', 0)
                        if 170 <= length <= 235:
                            seq = info.get('sequence', '')
                            alpha_candidates.append((cid, self._score_tcr_alpha(seq, length)))
                if alpha_candidates:
                    alpha_candidates.sort(key=lambda x: x[1], reverse=True)
                    best_cid, best_score = alpha_candidates[0]
                    if best_score >= 0.2:
                        assignments[best_cid] = 'tcr_alpha'
                    else:
                        # Length-based fallback: choose the shortest remaining unknown above 170 as alpha
                        shortest_unknown = min(
                            ((cid, info.get('length', 0)) for cid, info in chain_info.items() if assignments.get(cid, 'unknown') == 'unknown' and info.get('length', 0) >= 170),
                            key=lambda x: x[1],
                            default=None
                        )
                        if shortest_unknown:
                            assignments[shortest_unknown[0]] = 'tcr_alpha'
        except Exception:
            pass
        
        return assignments
    
    def _classify_multiple_complexes(self, complexes: List[Dict]) -> Dict[str, str]:
        """
        Classify chains within multiple MHC-II complexes.
        """
        all_assignments = {}
        
        for i, complex_chains in enumerate(complexes):
            if self.verbose:
                print(f"\nClassifying MHC-II Complex {i+1} ({len(complex_chains)} chains):")
                for chain_id, info in complex_chains.items():
                    print(f"  Chain {chain_id}: {len(info['sequence'])} residues")
            
            # Classify this complex
            complex_assignments = self._classify_chains(complex_chains)
            
            # Add complex number to chain assignments only if multiple complexes
            for chain_id, chain_type in complex_assignments.items():
                if len(complexes) > 1:
                    # Multiple complexes - add complex number
                    if chain_type != 'unknown':
                        all_assignments[chain_id] = f"{chain_type}_complex{i+1}"
                    else:
                        all_assignments[chain_id] = 'unknown'
                else:
                    # Single complex - use original names
                    all_assignments[chain_id] = chain_type
        
        return all_assignments
    
    # Scoring functions for each chain type
    
    def _score_peptide(self, sequence: str, length: int, properties: Dict) -> float:
        """Score likelihood of being a peptide antigen (MHC-II version - longer)."""
        score = 0.0
        
        # Length scoring for MHC-II peptides (longer than MHC-I)
        peptide_range = self.peptide_patterns['length_range']
        optimal_range = self.peptide_patterns['optimal_length']
        
        if peptide_range[0] <= length <= peptide_range[1]:
            score += 0.8  # Strong preference for correct length
            if optimal_range[0] <= length <= optimal_range[1]:
                score += 0.2  # Bonus for optimal length
        elif length < peptide_range[0]:
            # Too short for MHC-II peptide
            score -= 0.5
        elif length > peptide_range[1] + 5:
            # Too long for peptide
            score -= 0.8
        
        # Molecular weight check
        mw_range = self.peptide_patterns['molecular_weight']
        mw = properties.get('molecular_weight', 0)
        if mw_range[0] <= mw <= mw_range[1]:
            score += 0.2
        
        # Hydrophobic residues (common in peptide anchors)
        hydrophobic_count = sum(1 for aa in sequence if aa in self.peptide_patterns['common_residues'])
        hydrophobic_fraction = hydrophobic_count / length if length > 0 else 0
        score += hydrophobic_fraction * 0.3
        
        # Very short sequences are likely peptides
        if length <= 11:
            score += 0.4
        
        return min(max(score, 0.0), 1.0)
    
    def _score_mhc_ii_alpha(self, sequence: str, length: int) -> float:
        """Score likelihood of being MHC-II alpha chain."""
        score = 0.0
        
        # Length scoring
        alpha_range = self.mhc_ii_patterns['alpha']['length_range']
        optimal_range = self.mhc_ii_patterns['alpha']['optimal_length']
        
        if alpha_range[0] <= length <= alpha_range[1]:
            score += 0.7  # Strong weight for correct length
            if optimal_range[0] <= length <= optimal_range[1]:
                score += 0.3  # Bonus for optimal length
        elif length < alpha_range[0] - 10 or length > alpha_range[1] + 20:
            score -= 0.6  # Strong penalty for wrong length
        
        # Check for highly specific MHC-II α patterns
        n_terminal = sequence[:50] if len(sequence) >= 50 else sequence
        c_terminal = sequence[-50:] if len(sequence) >= 50 else sequence
        
        # Highly specific alpha motifs
        alpha_specific_found = 0
        for pattern in self.mhc_ii_patterns['alpha']['highly_specific_motifs']:
            if pattern in sequence:
                alpha_specific_found += 1
                score += 0.4  # High score for specific patterns
        
        # Multiple specific patterns = strong evidence
        if alpha_specific_found >= 2:
            score += 0.5
        
        # N-terminal signature patterns
        for pattern in self.mhc_ii_patterns['alpha']['n_terminal_signatures']:
            if pattern in n_terminal:
                score += 0.5  # Very strong MHC-II α indicator
                break
        
        # C-terminal patterns
        for pattern in self.mhc_ii_patterns['alpha']['c_terminal_patterns']:
            if pattern in c_terminal:
                score += 0.2
                break
        
        # Regular conserved motifs
        motifs = self.mhc_ii_patterns['alpha']['conserved_motifs']
        motif_count = sum(1 for motif in motifs if motif in sequence)
        score += motif_count * 0.1
        
        # Penalty for TCR-like patterns
        tcr_patterns = (self.tcr_patterns['alpha']['highly_specific_motifs'] + 
                       self.tcr_patterns['beta']['highly_specific_motifs'])
        if any(pattern in sequence for pattern in tcr_patterns):
            score -= 0.6
        
        return min(max(score, 0.0), 1.0)
    
    def _score_mhc_ii_beta(self, sequence: str, length: int) -> float:
        """Score likelihood of being MHC-II beta chain."""
        score = 0.0
        
        # Length scoring - more flexible for real sequences
        beta_range = self.mhc_ii_patterns['beta']['length_range']
        optimal_range = self.mhc_ii_patterns['beta']['optimal_length']
        
        if beta_range[0] <= length <= beta_range[1]:
            score += 0.6  # Strong weight for correct length
            if optimal_range[0] <= length <= optimal_range[1]:
                score += 0.2  # Bonus for optimal length
        elif 160 <= length <= 220:  # Extended acceptable range
            score += 0.4  # Partial credit for close length
        elif length < 150 or length > 250:
            score -= 0.5  # Penalty for clearly wrong length
        
        # Check for highly specific MHC-II β patterns
        n_terminal = sequence[:50] if len(sequence) >= 50 else sequence
        c_terminal = sequence[-50:] if len(sequence) >= 50 else sequence
        
        # Highly specific beta motifs
        beta_specific_found = 0
        for pattern in self.mhc_ii_patterns['beta']['highly_specific_motifs']:
            if pattern in sequence:
                beta_specific_found += 1
                score += 0.4  # High score for specific patterns
        
        # Multiple specific patterns = strong evidence
        if beta_specific_found >= 2:
            score += 0.5
        elif beta_specific_found >= 1:
            score += 0.3  # Some bonus for at least one pattern
        
        # N-terminal signature patterns
        for pattern in self.mhc_ii_patterns['beta']['n_terminal_signatures']:
            if pattern in n_terminal:
                score += 0.5  # Very strong MHC-II β indicator
                break
        
        # C-terminal patterns
        for pattern in self.mhc_ii_patterns['beta']['c_terminal_patterns']:
            if pattern in c_terminal:
                score += 0.2
                break
        
        # Regular conserved motifs
        motifs = self.mhc_ii_patterns['beta']['conserved_motifs']
        motif_count = sum(1 for motif in motifs if motif in sequence)
        score += motif_count * 0.1
        
        # General MHC-II patterns (less specific but more inclusive)
        general_mhc_patterns = [
            'FDSD',     # Common in MHC-II sequences
            'RVRL',     # Framework pattern
            'EYFM',     # Domain pattern
            'VYRA',     # Common motif
            'TRAELD',   # β2 domain pattern
            'CHVEH',    # β2 domain
        ]
        
        general_pattern_count = sum(1 for pattern in general_mhc_patterns if pattern in sequence)
        score += general_pattern_count * 0.15  # Moderate boost for general patterns
        
        # Amino acid composition scoring (MHC-II chains have characteristic composition)
        # High aromatic and charged residue content
        aromatic_count = sum(1 for aa in sequence if aa in 'FYW')
        charged_count = sum(1 for aa in sequence if aa in 'DEKR')
        aromatic_fraction = aromatic_count / length if length > 0 else 0
        charged_fraction = charged_count / length if length > 0 else 0
        
        if 0.05 <= aromatic_fraction <= 0.15:  # Typical aromatic content
            score += 0.1
        if 0.15 <= charged_fraction <= 0.35:  # Typical charged content
            score += 0.1
        
        # Fallback scoring for sequences with right length but no specific patterns
        if (160 <= length <= 220 and score < 0.3 and 
            beta_specific_found == 0):  # No specific patterns found
            # Check if it's clearly NOT another type
            not_peptide = length > 25
            not_tcr = all(all_scores.get('tcr_alpha', 0) < 0.5 and all_scores.get('tcr_beta', 0) < 0.5 
                         for all_scores in [{'tcr_alpha': self._score_tcr_alpha(sequence, length),
                                           'tcr_beta': self._score_tcr_beta(sequence, length)}])
            
            if not_peptide and not_tcr:
                score += 0.3  # Fallback score for process of elimination
        
        # Penalty for TCR-like patterns
        tcr_patterns = (self.tcr_patterns['alpha']['highly_specific_motifs'] + 
                       self.tcr_patterns['beta']['highly_specific_motifs'])
        tcr_pattern_count = sum(1 for pattern in tcr_patterns if pattern in sequence)
        if tcr_pattern_count > 0:
            score -= tcr_pattern_count * 0.3  # Scaled penalty
        
        return min(max(score, 0.0), 1.0)
    
    def _score_tcr_alpha(self, sequence: str, length: int) -> float:
        """Score likelihood of being TCR alpha chain (same as MHC-I analyzer)."""
        score = 0.0
        
        # Length check with stronger discrimination
        alpha_range = self.tcr_patterns['alpha']['length_range']
        optimal_range = self.tcr_patterns['alpha']['optimal_length']
        
        if alpha_range[0] <= length <= alpha_range[1]:
            score += 0.5  # Strong length preference
            # Bonus for optimal length range
            if optimal_range[0] <= length <= optimal_range[1]:
                score += 0.3
        elif 170 <= length < alpha_range[0]:
            # Slightly short alpha chains are acceptable in some PDBs
            score += 0.35
        elif length > alpha_range[1]:
            # Penalty for being too long (likely beta)
            score -= 0.4
        elif length < alpha_range[0]:
            # Penalty for being too short
            score -= 0.2
        
        # Check for highly specific alpha patterns
        n_terminal = sequence[:30] if len(sequence) >= 30 else sequence
        c_terminal = sequence[-30:] if len(sequence) >= 30 else sequence
        
        # Highly specific alpha motifs
        alpha_specific_found = 0
        for pattern in self.tcr_patterns['alpha']['highly_specific_motifs']:
            if pattern in sequence:
                alpha_specific_found += 1
                score += 0.4  # High score for specific patterns
        
        # Multiple specific patterns = strong evidence
        if alpha_specific_found >= 2:
            score += 0.4
        
        # N-terminal patterns (very distinctive for alpha)
        n_patterns = self.tcr_patterns['alpha']['n_terminal_patterns']
        for pattern in n_patterns:
            if pattern in n_terminal:
                score += 0.35
                break
        
        # C-terminal patterns
        c_patterns = self.tcr_patterns['alpha']['c_terminal_patterns']
        for pattern in c_patterns:
            if pattern in c_terminal:
                score += 0.25
                break
        
        # Regular motif checks
        motifs = self.tcr_patterns['alpha']['conserved_motifs']
        motif_count = sum(1 for motif in motifs if motif in sequence)
        score += motif_count * 0.1
        
        # Strong preference for shorter chains (alpha characteristic)
        if length < 190:
            score += 0.25
        elif length < 220:
            score += 0.2
        elif length > 250:
            score -= 0.5  # Strong penalty for long chains
        
        # Penalty for MHC-II-specific patterns
        mhc_alpha_patterns = self.mhc_ii_patterns['alpha']['highly_specific_motifs']
        mhc_beta_patterns = self.mhc_ii_patterns['beta']['highly_specific_motifs']
        if any(pattern in sequence for pattern in mhc_alpha_patterns + mhc_beta_patterns):
            score -= 0.8  # Strong penalty for MHC patterns
        
        # Penalty for beta-specific patterns
        beta_specific = self.tcr_patterns['beta']['highly_specific_motifs']
        if any(pattern in sequence for pattern in beta_specific):
            score -= 0.6  # Penalty for beta-specific patterns
        
        return min(max(score, 0.0), 1.0)
    
    def _score_tcr_beta(self, sequence: str, length: int) -> float:
        """Score likelihood of being TCR beta chain (same as MHC-I analyzer)."""
        score = 0.0
        
        # Length check with stronger discrimination
        beta_range = self.tcr_patterns['beta']['length_range']
        optimal_range = self.tcr_patterns['beta']['optimal_length']
        
        if beta_range[0] <= length <= beta_range[1]:
            score += 0.5  # Strong length preference
            # Bonus for optimal length range
            if optimal_range[0] <= length <= optimal_range[1]:
                score += 0.3
        elif 220 <= length < beta_range[0]:
            # Slightly short beta chains (common edge case)
            score += 0.2
        elif length < 220:
            # Penalty for being clearly too short (likely alpha)
            score -= 0.4
        elif length > beta_range[1]:
            # Penalty for being too long (likely MHC)
            score -= 0.3
        
        # Check for highly specific beta patterns
        n_terminal = sequence[:30] if len(sequence) >= 30 else sequence
        c_terminal = sequence[-30:] if len(sequence) >= 30 else sequence
        
        # Highly specific beta motifs
        beta_specific_found = 0
        for pattern in self.tcr_patterns['beta']['highly_specific_motifs']:
            if pattern in sequence:
                beta_specific_found += 1
                score += 0.4  # High score for specific patterns
        
        # Multiple specific patterns = strong evidence
        if beta_specific_found >= 2:
            score += 0.4
        
        # N-terminal patterns (very distinctive for beta)
        n_patterns = self.tcr_patterns['beta']['n_terminal_patterns']
        for pattern in n_patterns:
            if pattern in n_terminal:
                score += 0.35
                break
        
        # C-terminal patterns
        c_patterns = self.tcr_patterns['beta']['c_terminal_patterns']
        for pattern in c_patterns:
            if pattern in c_terminal:
                score += 0.25
                break
        
        # Regular motif checks
        motifs = self.tcr_patterns['beta']['conserved_motifs']
        motif_count = sum(1 for motif in motifs if motif in sequence)
        score += motif_count * 0.1
        
        # Strong preference for longer chains (beta characteristic)
        if length > 230 and length < 290:
            score += 0.2
        elif length < 220:
            score -= 0.5  # Strong penalty for short chains
        elif length > 300:
            score -= 0.6  # Very strong penalty for very long chains (likely MHC)
        
        # Penalty for MHC-II-specific patterns
        mhc_alpha_patterns = self.mhc_ii_patterns['alpha']['highly_specific_motifs']
        mhc_beta_patterns = self.mhc_ii_patterns['beta']['highly_specific_motifs']
        if any(pattern in sequence for pattern in mhc_alpha_patterns + mhc_beta_patterns):
            score -= 0.8  # Strong penalty for MHC patterns
        
        # Penalty for alpha-specific patterns
        alpha_specific = self.tcr_patterns['alpha']['highly_specific_motifs']
        if any(pattern in sequence for pattern in alpha_specific):
            score -= 0.6  # Penalty for alpha-specific patterns
        
        return min(max(score, 0.0), 1.0)
    
    def save_analysis_report(self, assignments: Dict[str, str], 
                           output_file: str, chain_info: Dict = None):
        """
        Save MHC-II analysis report to file.
        
        Args:
            assignments: Chain assignments
            output_file: Output file path
            chain_info: Additional chain information
        """
        with open(output_file, 'w') as f:
            f.write("MHC-II-TCR Complex Analysis Report\n")
            f.write("=" * 40 + "\n\n")
            
            # Sort by chain ID for consistent output
            for chain_id in sorted(assignments.keys()):
                chain_type = assignments[chain_id]
                f.write(f"Chain {chain_id}: {chain_type}\n")
                
                if chain_info and chain_id in chain_info:
                    info = chain_info[chain_id]
                    f.write(f"  Length: {info['length']} residues\n")
                    f.write(f"  Molecular Weight: {info['properties'].get('molecular_weight', 'N/A'):.1f} Da\n")
                    f.write(f"  Sequence: {info['sequence'][:50]}...\n\n")
        
        if self.verbose:
            print(f"Analysis report saved to: {output_file}")


def main():
    """Command-line interface for MHC-II analyzer."""
    import sys
    
    if len(sys.argv) != 2:
        print("Usage: python mhc_ii_analyzer.py <pdb_file>")
        sys.exit(1)
    
    pdb_file = sys.argv[1]
    
    try:
        analyzer = pMHCIITCRAnalyzer(verbose=True)
        assignments = analyzer.analyze_pdb(pdb_file)
        
        print("\nFinal Chain Assignments:")
        print("-" * 30)
        for chain_id in sorted(assignments.keys()):
            print(f"Chain {chain_id}: {assignments[chain_id]}")
        
    except Exception as e:
        print(f"Error: {e}")
        sys.exit(1)


if __name__ == "__main__":
    main()