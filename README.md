# TagGD: Barcode Demultiplexing Utilities for Spatial Transcriptomics Data

**TagGD** is a Python-based barcode demultiplexer for Spatial Transcriptomics data.
It provides a generalized, optimized, and up-to-date version of the original C++ demultiplexer "findIndexes," available [here](https://github.com/pelinakan/UBD).

For the original peer-reviewed reference to the program, see [PLOS ONE](http://www.plosone.org/article/info%3Adoi%2F10.1371%2Fjournal.pone.0057521).

## Overview

The primary goal of TagGD is to extract cDNA barcodes from input files (FASTQ, FASTA, SAM, or BAM)
and match them against a list of reference barcodes using a k-mer-based approach. Matched reads are
output with barcode and spatial information added to each record.

TagGD is versatile and can be used to demultiplex any type of index if a reference file is provided.
Users can even create fake spatial coordinates (X, Y) for general-purpose demultiplexing tasks.

## Key Features

- Supports FASTQ, FASTA, SAM, and BAM formats.
- Handles multiple indexes per read.
- K-mer-based matching for efficient and accurate demultiplexing.
- Outputs matched, unmatched, and ambiguous reads with annotated barcodes.
- Multiple options and distance metrice.
- Fast and memmory efficient.

---

## Requirements

- Python 3.10 or higher
- Cython
- PySam
- NumPy

---

## Installation

### From Source

If you are using a virtual environment like Anaconda:

```bash
git clone https://github.com/your-repo/taggd.git
cd taggd
python setup.py build
python setup.py install
```

### Using `pip`

Install directly from PyPI:

```bash
pip install taggd
```

---

## Building the Project

If you are contributing or making changes to the code, you may need to rebuild the Cython extensions:

```bash
python setup.py build_ext --inplace
```

## Testing the Project

```bash
pytest
```

---

## Usage

### Basic Command

To see all available options, run:

```bash
taggd_demultiplex.py -h
```

### Input Reference File Format

The reference file should contain barcodes and optional spatial coordinates, formatted as follows:

```tsv
BARCODE X Y
```

Example:

```tsv
ACGTACGT 0 0
TGCATGCA 1 1
```

---

### Example Commands

#### Example

```bash
taggd_demultiplex   --k 6   --max-edit-distance 3   --overhang 2   --subprocesses 4   --seed randomseed   <barcodes.tsv>   <input_file>   <output_prefix>
```

---

## Output

TagGD generates the following output files:

- `<output_prefix>_matched.*`: Reads that matched reference barcodes.
- `<output_prefix>_unmatched.*`: Reads that did not match any reference barcodes.
- `<output_prefix>_ambiguous.*`: Reads that matched multiple barcodes.
- `<output_prefix>_results.tsv`: Summary statistics of the run.

---

## Manual

### Options

Run `taggd_demultiplex -h` to view all available options and their descriptions.

---

## Contact

For questions, bug reports, or contributions, please contact:

- **Jose Fernandez Navarro**: <jc.fernandez.navarro@gmail.com>
