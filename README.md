# ARSC

[![PyPI version](https://badge.fury.io/py/arsc.svg)](https://badge.fury.io/py/arsc)
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.17745555.svg)](https://doi.org/10.5281/zenodo.17745555)
Visualization and retrieval of precomputed results can also be performed through [HERE](https://stsnsn.github.io/ARSC/).

**ARSC** (**A**mino-acid **R**esidue **S**toichiometry **C**alculator) is a lightweight command-line tool for quantifying elemental stoichiometry from protein FASTA files. It calculates the number of nitrogen (N), carbon (C), and sulfur (S) atoms contained in amino-acid side chains across all proteins, and also derives the average molecular weight of residues (AvgResMW).




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
#### 1. Run ARSC on a `.faa` file.
```bash
ARSC -i E_coli.faa
```

#### 2. Process all `.faa` / `.faa.gz` files in a directory using 4 threads and save results as `ARSC_output.tsv`
```bash
ARSC -i input_dir/ -t 4 -o ARSC_output.tsv
```

#### 3. Sort results by N-ARSC (descending) using pipe.
```bash
ARSC -i input_dir/ -t 4 | sort -k2,2nr
```

### Input requirements

- Input directory must contain one or more amino-acid sequence fasta (`*.faa` or `*.faa.gz`) files

### Output

- stdout
- TSV file (via `-o` or `--output`, optional)

Format columns: Genome, N_ARSC, C_ARSC, S_ARSC, AvgResMW <br>
- N-ARSC — Average number of nitrogen atoms per amino-acid residue side chain.
- C-ARSC — Average number of carbon atoms per amino-acid residue side chain.
- S-ARSC — Average number of sulfur atoms per amino-acid residue side chain.
- AvgResMW — Average molecular weight of amino-acid residues.

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