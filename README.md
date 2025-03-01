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
**Bioinformatics Utilities** is a **powerful and modular toolkit** for modern bioinformatics workflows.  
It includes utilities for **DNA/RNA sequence analysis, FASTQ filtering, and bioinformatics file processing**.

### Features  
✔️ **DNA/RNA Sequence Processing** – Transcription, complementarity, reversal, and GC content.  
✔️ **FASTQ Filtering** – Filter FASTQ sequences by GC content, length, and quality.  
✔️ **File Format Support** – FASTA, BLAST, and GBK file processing.  
✔️ **Optimized Performance** – Uses Biopython for high-speed computations.  

---

## Table of Contents

- [Installation](#installation)
- [Usage](#usage)
  - [1. Sequence Operations (run)](#1-sequence-operations-run)
  - [2. FASTQ Filtering (filter)](#2-fastq-filtering-filter)
  - [3. bio_files_processor.py](#3-bio_files_processorpy)
- [Project Structure](#project-structure)
- [Contributing](#contributing)

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

## <p align="center" style="color:#007acc;">Usage</p>

---

1. Sequence Operations (run)

The run subcommand in main_script.py enables you to perform various
operations on DNA/RNA sequences such as transcription, reversal, complementarity,
and GC content calculation.

Example Usage (Python API)

```python
from main_script import run_sequence_operation

sequences = ["ATGC", "CGTGA"]
operation = "reverse_complement"
result = run_sequence_operation("dna", sequences, operation)
print(result)
```

### <p align="center" style="color:#007acc;">Available Operations</p>

---

#### <p align="center" style="color:#009900;">DNA Sequences</p>

<p align="center">

| **Operation**            | **Description**                                      |
|-------------------------|------------------------------------------------------|
| `reverse`               | Reverses the sequence.                               |
| `complement`            | Returns the complementary sequence.                  |
| `reverse_complement`    | Returns the reverse complement.                      |
| `transcribe`            | Transcribes DNA to RNA.                              |

</p>

---

#### <p align="center" style="color:#ff6600;">RNA Sequences</p>

<p align="center">

| **Operation**            | **Description**                                      |
|-------------------------|------------------------------------------------------|
| `reverse`               | Reverses the sequence.                               |
| `complement`            | Returns the complementary sequence.                  |
| `reverse_complement`    | Returns the reverse complement.                      |

</p>

---

#### <p align="center" style="color:#cc0000;">Protein Sequences</p>

<p align="center">

| **Operation**            | **Description**                                      |
|-------------------------|------------------------------------------------------|
| `reverse`               | Reverses the sequence.                               |
| `amino_acid_composition` | Returns the amino acid composition.                 |
| `hydrophobicity_score`  | Calculates the hydrophobicity score.                 |

</p>

Command-Line Usage

```bash
python main_script.py run dna ATGC CGTGA reverse_complement
```

---

2. FASTQ Filtering (filter)

The filter subcommand in main_script.py filters FASTQ sequences
based on GC content, sequence length, and average quality.

Command-Line Usage

```bash
python main_script.py filter data/input.fastq filtered_output.fastq \
--gc_bounds 40 60 --length_bounds 50 1000 --quality_threshold 30
```

If the output_fastq argument is omitted, the tool returns a dictionary of
filtered records instead of writing to a file.

---

3. bio_files_processor.py

The bio_files_processor.py script provides utilities for converting and
processing bioinformatics file formats such as FASTA, BLAST, and GBK.

Key Functions

<p align="center">
  
| Command                          | Description                                      |
|-----------------------------------|--------------------------------------------------|
| **`convert_multiline_fasta_to_oneline`** | Converts multiline FASTA sequences to single-line. |
| **`parse_blast_output`**              | Extracts descriptions from BLAST output.        |
| **`select_genes_from_gbk_to_fasta`**  | Extracts genes from GBK and writes to FASTA.    |

</p>

#### Usage Examples:

#####  `convert_multiline_fasta_to_oneline`

```bash
python bio_files_processor.py convert_fasta input.fasta --output_fasta output.fasta
```
#####  parse_blast_output

```bash
python bio_files_processor.py parse_blast input_blast.txt output_descriptions.txt
```
#####  select_genes_from_gbk_to_fasta

```bash
python bio_files_processor.py select_genes input.gbk geneA geneB --n_before 1 --n_after 1 --output_fasta output.fasta
```

---

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

<p align="center">
  <img src="https://img.shields.io/badge/Quote-Einstein-blue?style=flat-square">
</p>

<p align="center"><em>“In the middle of difficulty lies opportunity.”</em></p>

<p align="center">— <strong>Albert Einstein</strong></p>

Enjoy using Bioinformatics Utilities and feel free to contribute to the project!
