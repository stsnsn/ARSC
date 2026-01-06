# -*- coding: utf-8 -*-
#!/usr/bin/env python3
__author__ = 'Satoshi_Nishino'
__email__ = 'satoshi-nishino@g.ecc.u-tokyo.ac.jp'


"""
This script was created to build a file processing utility for ARSC computations.
"""


import os
import re
import gzip
import subprocess
import shutil
import tempfile
import sys
from ARSC.core import process_faa  # compute_dna_metricsを削除
from collections import Counter
from Bio import SeqIO
from Bio.SeqUtils import gc_fraction


# Remove extensions
def get_genome_name(path):
    base = os.path.basename(path)
    root, ext = os.path.splitext(base)  # ext = ".faa" or ".gz"
    if ext == ".gz":
        root, _ = os.path.splitext(root)  # remove ".faa"
    return root


# Return items: {"handle": <file_path>, "name": <genome_name>}
def collect_faa_files(input_path):
    # --- directory ---
    if os.path.isdir(input_path):
        for f in os.listdir(input_path):
            if f.endswith(".faa") or f.endswith(".faa.gz"):
                fpath = os.path.join(input_path, f)
                genome = get_genome_name(f)
                yield {"handle": fpath, "name": genome}
        return

    # --- single file ---
    if os.path.isfile(input_path):
        if input_path.endswith(".faa") or input_path.endswith(".faa.gz"):
            genome =  get_genome_name(input_path)
            yield {"handle": input_path, "name": genome}
            return
        else:
            raise ValueError(f"input must be .faa or .faa.gz: {input_path}")

    raise ValueError(f"input path does not exist: {input_path}")


def collect_fna_files(input_path):
    # .fna, .fasta, .fa などを対象にする
    extensions = (".fna", ".fna.gz", ".fasta", ".fasta.gz", ".fa", ".fa.gz")
    if os.path.isdir(input_path):
        for f in os.listdir(input_path):
            if f.endswith(extensions):
                yield {"handle": os.path.join(input_path, f), "name": get_genome_name(f)}
    elif os.path.isfile(input_path) and input_path.endswith(extensions):
        yield {"handle": input_path, "name": get_genome_name(input_path)}


def process_faa_auto(item, per_sequence=False):
    """
    item: {"handle": path_str, "name": genome_name}
    per_sequence: bool, whether to process sequences individually
    """
    handle = item["handle"]
    name = item["name"]

    if handle.endswith(".gz"):
        with gzip.open(handle, "rt") as f:
            return process_faa(f, name=name, per_sequence=per_sequence)
    else:
        return process_faa(handle, name=name, per_sequence=per_sequence)


def run_prodigal(input_file):
    """
    Run Prodigal on the input file and return the path to the generated .faa file.
    """
    # Prodigalの出力ファイル名を定義
    name = get_genome_name(input_file)
    out_faa = os.path.join(tempfile.gettempdir(), f"{name}_output.faa")

    if input_file.endswith(".gz"):
        # 解凍してProdigal処理を行う
        with gzip.open(input_file, 'rb') as f_in:
            content = f_in.read()

        process = subprocess.Popen(
            ["prodigal", "-p", "meta", "-a", out_faa, "-o", "/dev/null"],
            stdin=subprocess.PIPE,
            stdout=subprocess.DEVNULL,
            stderr=subprocess.DEVNULL,
            text=True
        )
        process.communicate(input=content.decode())
    else:
        # 通常の.fnaファイル処理
        process = subprocess.Popen(
            ["prodigal", "-i", input_file, "-p", "meta", "-a", out_faa, "-o", "/dev/null"],
            stdout=subprocess.DEVNULL,
            stderr=subprocess.DEVNULL,
        )
        process.wait()

    return out_faa


def process_fna_pipeline(item, per_sequence=False):
    """
    Process .fna or .fna.gz files, handling GC content and Prodigal processing.
    """
    # Prodigalが存在しない場合のエラーハンドリング
    if shutil.which("prodigal") is None:
        print(f"Error: Prodigal not found in PATH. Skipping genome: {item['name']}", file=sys.stderr)
        return {"genome": item["name"], "error": "Prodigal not found in PATH."}

    handle = item["handle"]
    name = item["name"]

    base_counts = Counter()
    total_length = 0

    opener = gzip.open if handle.endswith(".gz") else open
    mode = "rt" if handle.endswith(".gz") else "r"
    
    try:
        with opener(handle, mode) as f:
            for line in f:
                if line.startswith(">"):
                    continue
                clean_line = line.strip().upper()
                base_counts.update(clean_line)
        
        # 集計
        A = base_counts.get("A", 0)
        T = base_counts.get("T", 0)
        G = base_counts.get("G", 0)
        C = base_counts.get("C", 0)
        total_atgc = A + T + G + C
        gc_content = (G + C) / total_atgc * 100 if total_atgc > 0 else 0
        
    except Exception as e:
        return {"genome": name, "error": f"Base composition error: {str(e)}"}

    # Prodigal
    faa_file = run_prodigal(handle)

    # ARSC 計算
    arsc_results = process_faa_auto({"handle": faa_file, "name": name}, per_sequence=per_sequence)

    # 結果
    result = {
        "GC": gc_content,
        "base_A": A*100/total_atgc if total_atgc > 0 else 0,
        "base_T": T*100/total_atgc if total_atgc > 0 else 0,
        "base_G": G*100/total_atgc if total_atgc > 0 else 0,
        "base_C": C*100/total_atgc if total_atgc > 0 else 0,
        **arsc_results
    }

    #(一時ファイル削除)
    if os.path.exists(faa_file):
        os.remove(faa_file)

    return result
