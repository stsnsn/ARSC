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
from multiprocessing import Pool
from ASTUR import __version__
from ASTUR.utils import collect_faa_files, process_faa_auto


ASTUR_LOGO = """
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;+;;xx+X;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;xX+xxxxxX++;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;xxxXXXXxXX+;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;+XXXXXXXxxX+;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;XXXXXXXXXx;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;xXXXXXXXXXx;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;xXXXXXXXXXX+;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;xXXXXXXXXXXx;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;+XXXXXXXXXXX+;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;xXXXXXXXXXXX+;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;xXXXXXXXXXXX+;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;+XXXXXXXXXXXx;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;+++;+XXXXXXXXXXXx;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;xXXXXXXXXXXXXXXXXx;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;xXXXXXXXXXXXX+;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;+XXXXXXXXXXXXXXXXXXXx;;;;;;;;
;;;;;;;;;;;;;;;;;;;+xXXXXXXXXXXXXXXXXXXXXx;;;;;;;;
;;;;;;;;;;;;;;;;;;xXXXXXXXXXXXXXXXXXXXXXXx;;;;;;;;
;;;;;;;;;;;;;;;;+XXXXXXXXXXXXXx+XXXXXXXXx+;;;;;;;;
;;;;;;;;;;;;;;+XXXXXXXXXXXXXXx;;xXXXXXXx+;;;;;;;;;
;;;;;;;;;;;;+XXXXXXXXXXXXXXX+;;;;XXXxX;;;;;;;;;;;;
;;;;;;;;;+xXXXXXXXXXXXXXXXx;;;;;;;+;;;;;;;;;;;;;;;
;;;;;;;+XXXXXXXXXXXXXXXXx;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;++XXXXXXXXXXXXx+;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;xxxxXXXXXXXx+;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;xxXxXxX++;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;xx+X++;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
 .d888888  .d88888b  d888888P dP     dP  888888ba 
d8'    88  88.    "'    88    88     88  88    `8b
88aaaaa88a `Y88888b.    88    88     88 a88aaaa8P'
88     88        `8b    88    88     88  88   `8b.
88     88  d8'   .8P    88    Y8.   .8P  88     88
88     88   Y88888P     dP    `Y88888P'  dP     dP
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
"""

class CustomFormatter(argparse.RawTextHelpFormatter, argparse.ArgumentDefaultsHelpFormatter):
    pass


def main():
    parser = argparse.ArgumentParser(description=f"{ASTUR_LOGO}\n\nCompute ARSC from .faa/.faa.gz files\nUsage example:\n    astur Ecoli.faa.gz\n    astur -i input_directory/ -o output.tsv -a -t 4 --stats", formatter_class=CustomFormatter)
    parser.add_argument("-i", "--input_dir", required=False, help="A faa or faa.gz file, or directory")
    parser.add_argument("input", nargs="?", help="Positional input: .faa/.faa.gz file or directory (or use -i/--input_dir)")
    parser.add_argument("-a", "--aa-composition", action="store_true", help="Include amino acid composition ratios and total length in output")
    parser.add_argument("-o", "--output", help="Output TSV file w/ header. If omitted, print to stdout w/ header.")
    parser.add_argument("--no-header", action="store_true", help="Suppress header line in output")
    parser.add_argument("-t", "--threads", default=1, type=int, help="Number of threads")
    parser.add_argument("-d", "--decimal-places", default=6, type=int, help="Number of decimal places for floating point values (default: 6)")
    parser.add_argument("--min-length", type=int, help="Minimum amino acid length (filter results)")
    parser.add_argument("--max-length", type=int, help="Maximum amino acid length (filter results)")
    parser.add_argument("-s", "--stats", action="store_true", help="Output summary statistics to stdout")
    parser.add_argument("-v", "--version", action="version", version=f"%(prog)s {__version__}")



    args = parser.parse_args()

    # Allow either flag (-i/--input_dir) or positional input (not both)
    if args.input_dir is not None and args.input is not None:
        parser.error("cannot specify both positional input and -i/--input_dir; use one or the other")
    if args.input_dir is None and args.input is not None:
        args.input_dir = args.input
    if args.input_dir is None:
        parser.error("missing input: provide a .faa/.faa.gz file or directory")

    # Check if input is a directory
    import os
    is_dir = os.path.isdir(args.input_dir)

    # --stats only works with directories
    if args.stats and not is_dir:
        parser.error("--stats can only be used with directory input (not single file)")

    try:
        items = list(collect_faa_files(args.input_dir))
    except ValueError as e:
        print(f"Error: {e}", file=sys.stderr)
        sys.exit(1)

    print(f"ASTUR Version: {__version__}", file=sys.stderr)
    print(f"Found {len(items)} files to process.", file=sys.stderr)
    print(f"Using {args.threads} threads.", file=sys.stderr)

    # Multiprocessing
    with Pool(args.threads) as pool:
        results = pool.map(process_faa_auto, items)

    # Filter by length
    filtered_results = []
    for r in results:
        if 'error' in r:
            filtered_results.append(r)
            continue
        length = r.get('total_aa_length', 0)
        if args.min_length is not None and length < args.min_length:
            continue
        if args.max_length is not None and length > args.max_length:
            continue
        filtered_results.append(r)
    
    results = filtered_results
    print(f"After filtering: {len(results)} results.", file=sys.stderr)

    # Calculate statistics if requested by --stats
    stats_data = None
    if args.stats:
        valid_results = [r for r in results if 'error' not in r]
        
        if valid_results:
            n_values = [r['N_ARSC'] for r in valid_results]
            c_values = [r['C_ARSC'] for r in valid_results]
            s_values = [r['S_ARSC'] for r in valid_results]
            mw_values = [r['MW_ARSC'] for r in valid_results]
            
            from statistics import mean, stdev
            
            stats_data = {
                'N_ARSC': {
                    'mean': mean(n_values),
                    'stdev': stdev(n_values) if len(n_values) > 1 else 0,
                    'min': min(n_values),
                    'max': max(n_values)
                },
                'C_ARSC': {
                    'mean': mean(c_values),
                    'stdev': stdev(c_values) if len(c_values) > 1 else 0,
                    'min': min(c_values),
                    'max': max(c_values)
                },
                'S_ARSC': {
                    'mean': mean(s_values),
                    'stdev': stdev(s_values) if len(s_values) > 1 else 0,
                    'min': min(s_values),
                    'max': max(s_values)
                },
                'MW_ARSC': {
                    'mean': mean(mw_values),
                    'stdev': stdev(mw_values) if len(mw_values) > 1 else 0,
                    'min': min(mw_values),
                    'max': max(mw_values)
                },
                'count': len(valid_results)
            }
            
            # Output statistics to stdout in a formatted table
            print("\n" + "="*70, file=sys.stderr)
            print("SUMMARY STATISTICS".center(70), file=sys.stderr)
            print("="*70, file=sys.stderr)
            print(f"{'Metric':<12} {'Mean':<16} {'Stdev':<16} {'Min':<16} {'Max':<16}", file=sys.stderr)
            print("-"*70, file=sys.stderr)
            for metric in ['N_ARSC', 'C_ARSC', 'S_ARSC', 'MW_ARSC']:
                s = stats_data[metric]
                print(f"{metric:<12} {s['mean']:<16.{args.decimal_places}f} {s['stdev']:<16.{args.decimal_places}f} {s['min']:<16.{args.decimal_places}f} {s['max']:<16.{args.decimal_places}f}", file=sys.stderr)
            print("-"*70, file=sys.stderr)
            print(f"{'Count':<12} {stats_data['count']:<16}", file=sys.stderr)
            print("="*70 + "\n", file=sys.stderr)

    # Output
    decimal_fmt = f"{{:.{args.decimal_places}f}}"

    if args.output:
        with open(args.output, "w") as out:
            if args.aa_composition:
                # Build header with all amino acids
                from ASTUR.core import aa_dictionary
                aa_keys = sorted(aa_dictionary.keys())
                # Write header unless --no-header is specified
                if not args.no_header:
                    header = "File\tN_ARSC\tC_ARSC\tS_ARSC\tAvgResMW\t" + "\t".join(aa_keys) + "\tTotalAALength\n"
                    out.write(header)
                for r in results:
                    if 'error' in r:
                        continue
                    aa_comp_values = [decimal_fmt.format(r['aa_composition'].get(aa, 0)) for aa in aa_keys]
                    out.write(f"{r['genome']}\t{decimal_fmt.format(r['N_ARSC'])}\t{decimal_fmt.format(r['C_ARSC'])}\t{decimal_fmt.format(r['S_ARSC'])}\t{decimal_fmt.format(r['MW_ARSC'])}\t" + "\t".join(aa_comp_values) + f"\t{r['total_aa_length']}\n")
            else:
                # Write header unless --no-header is specified
                if not args.no_header:
                    out.write("File\tN_ARSC\tC_ARSC\tS_ARSC\tAvgResMW\n")
                for r in results:
                    if 'error' in r:
                        continue
                    out.write(f"{r['genome']}\t{decimal_fmt.format(r['N_ARSC'])}\t{decimal_fmt.format(r['C_ARSC'])}\t{decimal_fmt.format(r['S_ARSC'])}\t{decimal_fmt.format(r['MW_ARSC'])}\n")
        print(f"Output written to {args.output}", file=sys.stderr)
    else:
        # stdout output with optional header
        if args.aa_composition:
            from ASTUR.core import aa_dictionary
            aa_keys = sorted(aa_dictionary.keys())
            # Print header unless --no-header is specified
            if not args.no_header:
                header = "File\tN_ARSC\tC_ARSC\tS_ARSC\tAvgResMW\t" + "\t".join(aa_keys) + "\tTotalAALength"
                print(header)
            for r in results:
                if 'error' in r:
                    continue
                aa_comp_values = [decimal_fmt.format(r['aa_composition'].get(aa, 0)) for aa in aa_keys]
                print(f"{r['genome']}\t{decimal_fmt.format(r['N_ARSC'])}\t{decimal_fmt.format(r['C_ARSC'])}\t{decimal_fmt.format(r['S_ARSC'])}\t{decimal_fmt.format(r['MW_ARSC'])}\t" + "\t".join(aa_comp_values) + f"\t{r['total_aa_length']}")
        else:
            # Print header unless --no-header is specified
            if not args.no_header:
                print("File\tN_ARSC\tC_ARSC\tS_ARSC\tAvgResMW")
            for r in results:
                if 'error' in r:
                    continue
                print(f"{r['genome']}\t{decimal_fmt.format(r['N_ARSC'])}\t{decimal_fmt.format(r['C_ARSC'])}\t{decimal_fmt.format(r['S_ARSC'])}\t{decimal_fmt.format(r['MW_ARSC'])}")



if __name__ == "__main__":
    main()