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

import sys
import argparse
import os
from multiprocessing import Pool
from statistics import mean, stdev
from ARSC import __version__
from ARSC.utils import collect_faa_files, process_faa_auto
from ARSC.core import aa_dictionary

quickARSC_LOGO = """
             _      _              _____   _____  _____ 
            (_)    | |       /\   |  __ \ / ____|/ ____|
  __ _ _   _ _  ___| | __   /  \  | |__) | (___ | |     
 / _` | | | | |/ __| |/ /  / /\ \ |  _  / \___ \| |     
| (_| | |_| | | (__|   <  / ____ \| | \ \ ____) | |____ 
 \__, |\__,_|_|\___|_|\_\/_/    \_\_|  \_\_____/ \_____|
    | |                                                 
    |_|                                                 
"""

class CustomFormatter(argparse.RawTextHelpFormatter, argparse.ArgumentDefaultsHelpFormatter):
    pass

def main():
    parser = argparse.ArgumentParser(description=f"{quickARSC_LOGO}\n\nCompute ARSC from .faa/.faa.gz files", formatter_class=CustomFormatter)
    parser.add_argument("input", help="Positional input: .faa/.faa.gz file or directory")
    parser.add_argument("-p", "--per-sequence", action="store_true", help="Process each sequence individually")
    parser.add_argument("-a", "--aa-composition", action="store_true", help="Include amino acid composition ratios")
    parser.add_argument("-o", "--output", help="Output TSV file (default: stdout)")
    parser.add_argument("--no-header", action="store_true", help="Suppress header line")
    parser.add_argument("-t", "--threads", default=1, type=int, help="Number of threads")
    parser.add_argument("-d", "--decimal-places", default=6, type=int, help="Decimal places")
    parser.add_argument("--min-length", type=int, help="Minimum sequence length")
    parser.add_argument("--max-length", type=int, help="Maximum sequence length")
    parser.add_argument("-s", "--stats", action="store_true", help="Output summary statistics to stderr")
    parser.add_argument("-v", "--version", action="version", version=f"%(prog)s {__version__}")

    args = parser.parse_args()

    # 入力パス
    target_input = args.input
    if not target_input:
        parser.error("missing input: provide a .faa/.faa.gz file or directory")

    if args.stats and not os.path.isdir(target_input):
        parser.error("--stats can only be used with directory input")

    try:
        items = list(collect_faa_files(target_input))
    except ValueError as e:
        print(f"Error: {e}", file=sys.stderr)
        sys.exit(1)

    print(f"ASTUR Version: {__version__}", file=sys.stderr)
    print(f"Found {len(items)} files to process.", file=sys.stderr)
    print(f"Using {args.threads} threads.", file=sys.stderr)

    # 計算
    with Pool(args.threads) as pool:
        if args.per_sequence:
            results = pool.starmap(process_faa_auto, [(item, True) for item in items])
        else:
            results = pool.map(process_faa_auto, items)

    # フィルタリング
    filtered_results = []
    for r in results:
        if 'error' in r:
            filtered_results.append(r)
            continue
        if args.per_sequence:
            r['sequences'] = [seq for seq in r['sequences'] if (args.min_length is None or seq['length'] >= args.min_length) and (args.max_length is None or seq['length'] <= args.max_length)]
            if r['sequences']:
                filtered_results.append(r)
        else:
            length = r.get('total_aa_length', 0)
            if (args.min_length is not None and length < args.min_length) or (args.max_length is not None and length > args.max_length):
                continue
            filtered_results.append(r)

    results = filtered_results
    print(f"After filtering: {len(results)} results.", file=sys.stderr)

    # 統計
    if args.stats:
        n_vals, c_vals, s_vals, mw_vals = [], [], [], []
        for r in results:
            if 'error' in r: continue
            target_list = r['sequences'] if args.per_sequence else [r]
            for data in target_list:
                n_vals.append(data.get('N_ARSC', 0))
                c_vals.append(data.get('C_ARSC', 0))
                s_vals.append(data.get('S_ARSC', 0))
                mw_vals.append(data.get('MW_ARSC', 0))
        
        if n_vals:
            print("\n" + "="*70, file=sys.stderr)
            label = "Per-Sequence" if args.per_sequence else "Per-File"
            print(f"SUMMARY STATISTICS ({label})".center(70), file=sys.stderr)
            print("="*70, file=sys.stderr)
            print(f"{'Metric':<12} {'Mean':<16} {'Stdev':<16} {'Min':<16} {'Max':<16}", file=sys.stderr)
            print("-"*70, file=sys.stderr)
            for name, vals in [('N_ARSC', n_vals), ('C_ARSC', c_vals), ('S_ARSC', s_vals), ('AvgResMW', mw_vals)]:
                sd = stdev(vals) if len(vals) > 1 else 0
                print(f"{name:<12} {mean(vals):<16.{args.decimal_places}f} {sd:<16.{args.decimal_places}f} {min(vals):<16.{args.decimal_places}f} {max(vals):<16.{args.decimal_places}f}", file=sys.stderr)
            print("-"*70, file=sys.stderr)
            print(f"{'Count':<12} {len(n_vals):<16}", file=sys.stderr)
            print("="*70 + "\n", file=sys.stderr)

    # 出力
    decimal_fmt = f"{{:.{args.decimal_places}f}}"
    out_handle = open(args.output, "w") if args.output else sys.stdout
    aa_keys = sorted(aa_dictionary.keys())

    try:
        # ヘッダー
        if not args.no_header:
            h = ["query"]
            if args.per_sequence:
                h.append("sequence_id")
            h.extend(["N_ARSC", "C_ARSC", "S_ARSC", "AvgResMW", "length" if args.per_sequence else "TotalLength"])
            if args.aa_composition:
                h.extend(aa_keys)
            out_handle.write("\t".join(h) + "\n")

        # データ行
        for r in results:
            if 'error' in r: continue
            
            if args.per_sequence:
                for seq in r.get('sequences', []):
                    row = [r['genome'], seq['sequence_id']]
                    row.extend([
                        decimal_fmt.format(seq.get('N_ARSC', 0)),
                        decimal_fmt.format(seq.get('C_ARSC', 0)),
                        decimal_fmt.format(seq.get('S_ARSC', 0)),
                        decimal_fmt.format(seq.get('MW_ARSC', 0)),
                        str(seq.get('length', 0))
                    ])
                    if args.aa_composition:
                        comp = seq.get('aa_composition', {})
                        row.extend([decimal_fmt.format(comp.get(aa, 0)) for aa in aa_keys])
                    out_handle.write("\t".join(row) + "\n")
            else:
                row = [r['genome']]
                row.extend([
                    decimal_fmt.format(r.get('N_ARSC', 0)),
                    decimal_fmt.format(r.get('C_ARSC', 0)),
                    decimal_fmt.format(r.get('S_ARSC', 0)),
                    decimal_fmt.format(r.get('MW_ARSC', 0)),
                    str(r.get('total_aa_length', 0))
                ])
                if args.aa_composition:
                    comp = r.get('aa_composition', {})
                    row.extend([decimal_fmt.format(comp.get(aa, 0)) for aa in aa_keys])
                out_handle.write("\t".join(row) + "\n")

    finally:
        if args.output:
            out_handle.close()
            print(f"Output written to {args.output}", file=sys.stderr)

if __name__ == "__main__":
    main()