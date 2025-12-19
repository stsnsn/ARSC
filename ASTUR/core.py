# -*- coding: utf-8 -*-
#!/usr/bin/env python3
__author__ = 'Satoshi_Nishino'
__email__ = 'satoshi-nishino@g.ecc.u-tokyo.ac.jp'


"""
This script was created to compute the N/C-ARSC metrics as described in the following publication:
    Mende et al., Nature Microbiology, 2017 https://doi.org/10.1038/s41564-017-0008-3

Original citations for calculation metrics in Mende et al. 2017:
    Baudouin-Cornu P, Surdin-Kerjan Y, Marliere P, Thomas D. 2001. Molecular evolution of protein atomic composition. Science 293 297–300.
    Wright F. 1990. The 'effective number of codons' used in a gene. Gene 87 23-29.
"""


import os
import sys
from Bio import SeqIO
from collections import Counter

# Amino acid dictionary
# -----------------------
aa_dictionary = {
    'K': {'N': 1, 'S': 0, 'MW': 146.1882, 'C': 4},
    'R': {'N': 3, 'S': 0, 'MW': 174.2017, 'C': 4},
    'H': {'N': 2, 'S': 0, 'MW': 155.1552, 'C': 4},
    'D': {'N': 0, 'S': 0, 'MW': 133.1032, 'C': 2},
    'E': {'N': 0, 'S': 0, 'MW': 147.1299, 'C': 3},
    'N': {'N': 1, 'S': 0, 'MW': 132.1184, 'C': 2},
    'Q': {'N': 1, 'S': 0, 'MW': 146.1451, 'C': 3},
    'S': {'N': 0, 'S': 0, 'MW': 105.0930, 'C': 1},
    'T': {'N': 0, 'S': 0, 'MW': 119.1197, 'C': 2},
    'Y': {'N': 0, 'S': 0, 'MW': 181.1894, 'C': 7},
    'A': {'N': 0, 'S': 0, 'MW': 89.0935,  'C': 1},
    'V': {'N': 0, 'S': 0, 'MW': 117.1469, 'C': 3},
    'L': {'N': 0, 'S': 0, 'MW': 131.1736, 'C': 4},
    'I': {'N': 0, 'S': 0, 'MW': 131.1736, 'C': 4},
    'P': {'N': 0, 'S': 0, 'MW': 115.1310, 'C': 3},
    'F': {'N': 0, 'S': 0, 'MW': 165.1900, 'C': 7},
    'M': {'N': 0, 'S': 1, 'MW': 149.2124, 'C': 3},
    'W': {'N': 1, 'S': 0, 'MW': 204.2262, 'C': 9},
    'G': {'N': 0, 'S': 0, 'MW': 75.0669,  'C': 0},
    'C': {'N': 0, 'S': 1, 'MW': 121.1590, 'C': 1},
    'U': {'N': 0, 'S': 0, 'MW': 168.07,   'C': 1},
    'J': {'N': 0, 'S': 0, 'MW': 131.1736, 'C': 4},
    'B': {'N': 0.5, 'S': 0, 'MW': 132.6108, 'C': 2},
    'Z': {'N': 0.5, 'S': 0, 'MW': 146.6375, 'C': 3}
}

# Compute ARSCs (N/C/S/MW)
# -----------------------
def compute_ARSC_extended_counts(counts, aa_dict):
    total_aa = sum(counts.values())
    if total_aa == 0:
        return None, None, None, None

    # Detect and log ignored characters
    ignored_chars = [a for a in counts if a not in aa_dict]
    if ignored_chars:
        print(f"Warning: Ignored characters found in sequence: {', '.join(set(ignored_chars))}", file=sys.stderr)

    total_N  = sum(counts[a] * aa_dict[a]["N"]  for a in counts if a in aa_dict)
    total_C  = sum(counts[a] * aa_dict[a]["C"]  for a in counts if a in aa_dict)
    total_S  = sum(counts[a] * aa_dict[a]["S"]  for a in counts if a in aa_dict)
    total_MW = sum(counts[a] * aa_dict[a]["MW"] for a in counts if a in aa_dict)

    return (
        total_N  / total_aa,
        total_C  / total_aa,
        total_S  / total_aa,
        total_MW / total_aa
    )


def compute_aa_composition(counts):
    """Compute amino acid composition ratios."""
    total_aa = sum(counts.values())
    if total_aa == 0:
        return {}
    
    # Return composition for all amino acids in aa_dictionary
    composition = {}
    for aa in aa_dictionary.keys():
        composition[aa] = counts.get(aa, 0) / total_aa
    return composition


def process_faa(faa_source, name=None, per_sequence=False):
    try:
        # Determine genome name
        genome_name = name
        if genome_name is None:
            # faa_source is a path string
            base = os.path.basename(faa_source)
            for ext in [".faa.gz", ".faa", ".gz"]:
                if base.endswith(ext):
                    base = base[: -len(ext)]
            genome_name = base

        if per_sequence:
            results = []
            for record in SeqIO.parse(faa_source, "fasta-blast"):
                seq = str(record.seq).upper().replace("*", "")
                seq_counts = Counter(seq)
                seq_length = sum(seq_counts.values())
                N, C, S, MW = compute_ARSC_extended_counts(seq_counts, aa_dictionary)
                results.append({
                    "sequence_id": record.id,
                    "length": seq_length,
                    "N_ARSC": N,
                    "C_ARSC": C,
                    "S_ARSC": S,
                    "MW_ARSC": MW,
                    "aa_composition": compute_aa_composition(seq_counts)
                })
            return {"genome": genome_name, "sequences": results}
        else:
            counts = Counter()
            for record in SeqIO.parse(faa_source, "fasta-blast"):
                seq = str(record.seq).upper().replace("*", "")
                counts.update(seq)

            total_aa_length = sum(counts.values())
            N, C, S, MW = compute_ARSC_extended_counts(counts, aa_dictionary)
            aa_composition = compute_aa_composition(counts)

            return {
                "genome": genome_name,
                "N_ARSC": N,
                "C_ARSC": C,
                "S_ARSC": S,
                "MW_ARSC": MW,
                "aa_composition": aa_composition,
                "total_aa_length": total_aa_length
            }

    except Exception as e:
        return {"genome": genome_name if name else None, "error": str(e)}
