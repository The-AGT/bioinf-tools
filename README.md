# Bioinformatics Utilities

[![Build Status](https://img.shields.io/badge/build-passing-brightgreen.svg)](https://github.com/The-AGT/bioinf-utils)
[![Flake8 Compliance](https://img.shields.io/badge/flake8-compliant-brightgreen.svg)](https://flake8.pycqa.org)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](LICENSE)
[![Python Version](https://img.shields.io/badge/Python-3.8%2B-blue.svg)](https://www.python.org/)

---

<p align="center">
  <img src="https://raw.githubusercontent.com/The-AGT/bioinf-tools/bioinf_tools_v2/utils/screenshots/cover_image.png" alt="Bioinformatics Utilities" width="600">
</p>

---

## Overview

**Bioinformatics Utilities** is a comprehensive toolkit designed for modern  
bioinformatics workflows. It provides utilities for:

- **DNA/RNA Sequence Processing:**  
  Transcription, complementarity, reverse, and GC content calculation.
- **FASTQ Filtering:**  
  Filter FASTQ files based on GC content, sequence length, and quality.
- **Bioinformatics File Processing:**  
  Convert and process FASTA, BLAST, and GBK file formats.

Leveraging industry-standard libraries such as Biopython, this toolkit  
ensures reliability and performance for research and production environments.

---

## Table of Contents

- [Installation](#installation)
- [Usage](#usage)
  - [1. run_dna_rna_tools](#1-run_dna_rna_tools)
  - [2. FASTQ Filtering](#2-fastq-filtering)
  - [3. bio_files_processor.py](#3-bio_files_processorpy)
  - [4. filter_fastq (legacy)](#4-filter_fastq-legacy)
- [Project Structure](#project-structure)
- [Screenshots](#screenshots)
- [Contributing](#contributing)
- [License](#license)

---

## Installation

### Clone the Repository

```bash
git clone https://github.com/The-AGT/bioinf-utils.git

Navigate to the Project Directory

cd bioinf-utils

Create a Virtual Environment (Recommended)

python -m venv venv
source venv/bin/activate   # Windows: venv\Scripts\activate

Install Dependencies

pip install -r requirements.txt

Usage

1. run_dna_rna_tools

The run_dna_rna_tools function (in main_script.py) allows you to
perform operations on DNA/RNA sequences such as transcription, reversal,
complementarity, and GC content calculation.

Example Usage

from main_script import run_dna_rna_tools

sequences = ["ATGC", "CGTGA"]
procedure = "reverse_complement"
result = run_dna_rna_tools(*sequences, procedure)
print(result)

Available Procedures
	•	transcribe
Transcribes DNA into RNA (T → U).
	•	reverse
Reverses the sequence.
	•	complement
Returns the complementary sequence.
	•	reverse_complement
Returns the reverse complement.
	•	gc_content
Calculates the GC content.

2. FASTQ Filtering

FASTQ filtering is now fully integrated into main_script.py via the
filter subcommand. This tool filters FASTQ sequences based on GC content,
sequence length, and quality.

Command-Line Usage

python main_script.py filter data/input.fastq filtered_output.fastq \
--gc_bounds 40 60 --length_bounds 50 1000 --quality_threshold 30

If the output_fastq argument is omitted, the tool returns a dictionary
of filtered records instead of writing to a file.

3. bio_files_processor.py

The bio_files_processor.py script offers utilities for converting and
processing bioinformatics file formats such as FASTA, BLAST, and GBK.

Key Functions
	1.	convert_multiline_fasta_to_oneline
Converts FASTA files with multiline sequences to a one-line format.
Usage:

python bio_files_processor.py convert_fasta input.fasta --output_fasta output.fasta


	2.	parse_blast_output
Extracts descriptions of the best matches from a BLAST output file.
Usage:

python bio_files_processor.py parse_blast input_blast.txt output_descriptions.txt


	3.	select_genes_from_gbk_to_fasta
Extracts protein sequences from a GBK file for target genes and their
neighbors, and writes them to a FASTA file.
Usage:

python bio_files_processor.py select_genes input.gbk geneA geneB \
--n_before 1 --n_after 1 --output_fasta output.fasta

4. filter_fastq (Legacy)

The legacy filter_fastq function supports filtering on a dictionary of
FASTQ records (key: name, value: tuple(sequence, quality, extra_info)).

Example Usage

from main_script import filter_fastq

seqs = {
    "seq1": ("ATGC", "IIII", "extra_info1"),
    "seq2": ("CGTGA", "HHHHH", "extra_info2")
}
gc_bounds = (40, 60)
length_bounds = (4, 10)
quality_threshold = 30

filtered_seqs = filter_fastq(seqs, gc_bounds, length_bounds, quality_threshold)
print(filtered_seqs)

Project Structure

bioinf-utils/
├── README.md
├── main_script.py
├── bio_files_processor.py
├── requirements.txt
└── utils/
    ├── __init__.py
    ├── dna_rna_tools_test.py
    ├── fastq_filtrator_test.py
    └── example_data.py

Screenshots

<p align="center">
  <img src="https://github.com/The-AGT/bioinf-tools/blob/bioinf_tools_v2/utils/screenshots/HW4_flake8.png" 
       alt="Flake8 Test Results" width="600">
</p>


Contributing

Contributions are welcome! Please fork this repository and submit
a pull request for any improvements or bug fixes.

For major changes, please open an issue first to discuss what you would like to change.

License

This project is licensed under the MIT License.

	“Innovation distinguishes between a leader and a follower.”
— Steve Jobs