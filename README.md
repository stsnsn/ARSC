# quickARSC: ARSC-based stoichiometry utility

[![PyPI version](https://badge.fury.io/py/arsc.svg)](https://badge.fury.io/py/arsc)
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.17977582.svg)](https://doi.org/10.5281/zenodo.17977582)


---

**quickARSC** is [a lightweight command-line tool](https://pypi.org/project/arsc/) and [a web interface](https://stsnsn.github.io/quickARSC/) for quantifying elemental stoichiometry from protein FASTA files. It calculates the number of nitrogen (N), carbon (C), and sulfur (S) atoms per amino acid residue side chain (ARSC) across proteins or proteomes.

These metrics follow the definitions used in
Mende et al., *Nature Microbiology*, (2017). https://doi.org/10.1038/s41564-017-0008-3

---

## Web Interface

The static web interface is available at [https://stsnsn.github.io/quickARSC/](https://stsnsn.github.io/quickARSC/)

**Features:**
- **Pre-computed Results**: Browse and download ARSC metrics for all 143,614 GTDB r226.0 representatives.
- **Interactive Filtering**: Filter results by taxonomy information.
- **Custom Analysis**: Upload your own amino acid FASTA file (.fa, .faa, .fasta) to compute ARSC metrics on-the-fly.

---

## Standalone Package

The standalone package is available at [https://pypi.org/project/arsc/](https://pypi.org/project/arsc/)

**Features**
- **Elemental stoichiometry calculation**: Calculate N-, C-, and S-ARSC directly from protein FASTA files or directories.
- **Multiprocessing**: Fast and scalable analysis of large genome or proteome datasets.
- **Simple CLI tool**: Run with a single command; easy to integrate into UNIX pipelines.

---

### Installation From PyPI

```bash
pip install arsc
```

---

## Usage

```bash
arsc <FASTA_FILE (.faa/.faa.gz) or input_dir/>
```

- `-h` or `--help`    : show help message
- `-v` or `--version` : show version

- `-o` or `--output`   <output> : output TSV file name (optional)
- `-t` or `--threads` N : number of threads (default: 1)
- `-s` or `--stats`     : output summary statistics to stderr (default: False)
- `-p`, `--per-sequence`: process each sequence individually instead of the entire file

- output format options
    - `-a` or `--aa-composition`   : Include amino acid composition ratios in output (default: False)
    - `-d` or `--decimal-places` N : Number of decimal places (default: 6)
    - `--no-header`    : Suppress header line in output (default: False)
    - `--max-length` N : number of maximal amino acid length (default: None)
    - `--min-length` N : number of minimal amino acid length (default: None)


### Example
#### 1. Compute ARSC values on a `.faa` file.
```bash
arsc test_data/genome_a.faa
```
- output example:

| query | N_ARSC | C_ARSC | S_ARSC | AvgResMW | TotalLength |
| --- | --- | --- | --- | --- | --- |
| genome_a | 0.148438 | 3.132812 | 0.023438 | 123.568566 | 194 |

#### 2. Process all `.faa` / `.faa.gz` files in a directory using 3 threads and save results as `ARSC_output.tsv`.
```bash
arsc test_data/ -t 3 -o ARSC_output.tsv
```

#### 3. Output with amino acid composition table as `ARSC_output_full.tsv` and show statistics summary.
```bash
arsc test_data/ -t 3 -as -o ARSC_output_full.tsv
```

#### 4. Sort results by N-ARSC (descending) using pipe.
```bash
arsc test_data/ -t 3 --no-header | sort -k2,2nr
```

#### 5. Process each sequence individually instead of the entire file and filter results by amino acid length > 120.
```bash
arsc test_data/ -t 3 --min-length 120 -p
```

### Input requirements

- Input directory must contain one or more amino-acid fasta (`*.faa` or `*.faa.gz`) files

### Output

- stdout (if you need no header, use `--no-header` option)
- TSV file (via `-o` or `--output`, optional)

Default format columns: query, N_ARSC, C_ARSC, S_ARSC, AvgResMW, TotalLenghth <br>
- N-ARSC — Average number of nitrogen atoms per amino-acid residue side chain.
- C-ARSC — Average number of carbon atoms per amino-acid residue side chain.
- S-ARSC — Average number of sulfur atoms per amino-acid residue side chain.
- AvgResMW — Average molecular weight of amino-acid residues (not only side chain!).
- TotalLenghth — Total amino acid length.

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
