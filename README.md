# ARSC

[![PyPI version](https://badge.fury.io/py/arsc.svg)](https://badge.fury.io/py/arsc)

ARSC (Amino-acid Residue Stoichiometric Constraints) is a lightweight tool to compute
**N-ARSC**, **C-ARSC**, **S-ARSC**, and **AvgResMW** from genomic protein FASTA files.

These metrics follow the definitions used in
Mende et al., *Nature Microbiology*, (2017). https://doi.org/10.1038/s41564-017-0008-3


---

## Features

- Calculate elemental composition metrics (N-ARSC, C-ARSC, S-ARSC, AvgResMW) directly from protein FASTA files.
- Multiprocessing support for fast and scalable analysis of large genome sets.
- Simple CLI tool: one command to run, easy to combine with UNIX tools via pipes.

---

## Installation

### From PyPI

```bash
pip install arsc
```

---

## Usage

```bash
ARSC -i <input> -o <output.tsv> -t <num_threads>
```

- `-v` or `--version` : show version
- `-h` or `--help`    : show help message

- `-i` : input file/directory path
- `-o` : output TSV file name (optional)
- `-t` : number of threads (default: 1)


### Example
#### 1. Calculate ARSC from a `.faa` file and save result as `ARSC_output.tsv`.
```bash
ARSC -i E_coli.faa -o ARSC_output.tsv
```

#### 2. Calculate ARSC for all .faa / .faa.gz files in a directory using 4 threads
```bash
ARSC -i input_dir/ -t 4
```

#### 3. Sort the result by N-ARSC.
```bash
ARSC -i input_dir/ -t 4 | sort -k2,2nr
```

### Input requirements

- Input directory must contain one or more amino-acid sequence fasta (`*.faa` or `*.faa.gz`) files

### Output format
- stdout
- TSV (-o or --output, optional)
    - containing Genome, N_ARSC, C_ARSC, S_ARSC, and AvgResMW.

### Dependencies
- Python >= 3.8
- Biopython >= 1.79

---

## Citation
Please cite following articles:
- (To be added)
- Mende et al., *Nature Microbiology*, (2017). https://doi.org/10.1038/s41564-017-0008-3

---

## License
This project is distributed under the GPL-2.0 license.

---