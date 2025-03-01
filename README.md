# <span style="color:#0056b3;">Bioinformatics Utilities</span>

[![Build Status](https://img.shields.io/badge/build-passing-brightgreen.svg)](https://github.com/The-AGT/bioinf-utils)
[![Flake8 Compliance](https://img.shields.io/badge/flake8-compliant-brightgreen.svg)](https://flake8.pycqa.org)
[![Python Version](https://img.shields.io/badge/Python-3.8%2B-blue.svg)](https://www.python.org/)

---
<p align="center">
  <img src="https://img.shields.io/badge/Bioinformatics-Utilities-blue?style=for-the-badge&logo=dna" alt="Bioinformatics Utilities">
</p>
---

## Overview

**Bioinformatics Utilities** is a comprehensive toolkit designed for modern  
bioinformatics workflows. It provides utilities for:

- **DNA/RNA Sequence Processing:**  
  Transcription, complementarity, reversal, and GC content calculation.
- **FASTQ Filtering:**  
  Filter FASTQ files based on GC content, sequence length, and quality.
- **Bioinformatics File Processing:**  
  Conversion and processing of FASTA, BLAST, and GBK file formats.

Leveraging industry-standard libraries such as Biopython, this toolkit  
ensures reliability and performance for both research and production environments.

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
```

Navigate to the Project Directory

```bash
cd bioinf-utils
```

Create a Virtual Environment (Recommended)

```bash
python -m venv venv
source venv/bin/activate   # For Windows: venv\Scripts\activate
```

Install Dependencies

```bash
pip install -r requirements.txt
```

Usage

1. run_dna_rna_tools

The run_dna_rna_tools function (located in main_script.py)
allows you to perform operations on DNA/RNA sequences such as transcription,
reversal, complementarity, and GC content calculation.

Example Usage

```python
from main_script import run_dna_rna_tools

sequences = ["ATGC", "CGTGA"]
procedure = "reverse_complement"
result = run_dna_rna_tools(*sequences, procedure)
print(result)
```

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

FASTQ filtering functionality is integrated into main_script.py via the
filter subcommand. This tool filters FASTQ sequences based on GC content,
sequence length, and quality.

Command-Line Usage

```bash
python main_script.py filter data/input.fastq filtered_output.fastq \
--gc_bounds 40 60 --length_bounds 50 1000 --quality_threshold 30
```

If the output_fastq argument is omitted, the tool returns a dictionary of
filtered records instead of writing to a file.

3. bio_files_processor.py

The bio_files_processor.py script provides utilities for converting and
processing bioinformatics file formats such as FASTA, BLAST, and GBK.

Key Functions

	a.	convert_multiline_fasta_to_oneline
Converts a FASTA file with multiline sequences into one with single-line
sequences.
Usage:

```bash
python bio_files_processor.py convert_fasta input.fasta --output_fasta output.fasta
```

	b.	parse_blast_output
Parses a BLAST output file and extracts descriptions of the best matches.
Usage:

```bash
python bio_files_processor.py parse_blast input_blast.txt output_descriptions.txt
```

	c.	select_genes_from_gbk_to_fasta
Extracts protein sequences from a GBK file for specified genes and their
neighbors, and writes them to a FASTA file.
Usage:

```bash
python bio_files_processor.py select_genes input.gbk geneA geneB \
--n_before 1 --n_after 1 --output_fasta output.fasta
```

4. filter_fastq (Legacy)

For users preferring a dictionary of sequences over a FASTQ file, the
legacy filter_fastq function is available. It operates on a dictionary
where each key is a sequence name and the value is a tuple
(sequence, quality, extra_info).

Example Usage

```python
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
```

Project Structure
```bash
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
```

Contributing

Contributions are welcome! Please fork this repository and submit a pull
request for any improvements or bug fixes.

For major changes, please open an issue first to discuss what you would like
to change.

	“In the middle of difficulty lies opportunity.”
— Albert Einstein

Enjoy using Bioinformatics Utilities and feel free to contribute to the project!
