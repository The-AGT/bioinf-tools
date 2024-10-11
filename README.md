```markdown
# Bioinformatics Utilities

## Description

This project provides a set of utilities for working with DNA and RNA sequences, filtering FASTQ files, and processing various bioinformatics data formats such as FASTA, BLAST, and GBK files. The main functionalities of the project include:

- **run_dna_rna_tools** — a set of tools for working with DNA/RNA sequences (transcription, complementarity, reverse, and other operations).
- **filter_fastq_file** — a utility for filtering FASTQ sequences based on GC content, sequence length, and quality. Functions for reading and writing FASTQ files are moved into a separate module (`utils/fastq_module.py`).
- **bio_files_processor.py** — a set of utilities for converting and processing bioinformatics file formats (FASTA, BLAST, and GBK).

## Installation

1. Clone the repository:
    ```bash
    git clone https://github.com/The-AGT/bioinf-utils.git
    ```

2. Navigate to the project directory:
    ```bash
    cd bioinf-utils
    ```

3. Create a virtual environment (recommended):
    ```bash
    python -m venv venv
    source venv/bin/activate  # For Windows: venv\Scripts\activate
    ```

4. Install dependencies (if any):
    ```bash
    pip install -r requirements.txt
    ```

## Usage

### 1. run_dna_rna_tools

The `run_dna_rna_tools` function, located in `main_script.py`, allows you to perform operations on DNA and RNA sequences, such as transcription, reverse, complementarity, and GC content calculation.

#### Example usage:
```python
from main_script import run_dna_rna_tools

sequences = ["ATGC", "CGTGA"]
procedure = "reverse_complement"
result = run_dna_rna_tools(*sequences, procedure)
print(result)
```

#### Available procedures:
- `transcribe` — transcribes DNA into RNA (replaces T with U).
- `reverse` — reverses the sequence.
- `complement` — returns the complementary sequence.
- `reverse_complement` — returns the reverse complement of the sequence.
- `gc_content` — calculates the GC content of the sequence.

### 2. filter_fastq_file

The `filter_fastq_file` function is used to filter FASTQ sequences based on GC content, sequence length, and average quality. This function works directly with FASTQ files, reading them, filtering the sequences on-the-fly, and writing the filtered sequences to an output FASTQ file. The functions for reading and writing FASTQ files have been moved into `utils/fastq_module.py`.

#### Arguments:
- **input_fastq** — path to the input FASTQ file.
- **output_fastq** — name of the output FASTQ file.
- **gc_bounds** — GC content bounds for filtering (default is (0, 100)).
- **length_bounds** — sequence length bounds for filtering (default is (0, 2**32)).
- **quality_threshold** — threshold for the average quality score (default is 0).

#### Example usage:
```python
from utils.fastq_module import filter_fastq_file

input_fastq = "data/input.fastq"
output_fastq = "filtered_output.fastq"
gc_bounds = (40, 60)
length_bounds = (50, 1000)
quality_threshold = 30

filter_fastq_file(input_fastq, output_fastq, gc_bounds, length_bounds, quality_threshold)
```

### 3. bio_files_processor.py

The `bio_files_processor.py` script provides utilities for converting and processing bioinformatics file formats, including FASTA, BLAST, and GBK files.

#### Functions in `bio_files_processor.py`:

1. **convert_multiline_fasta_to_oneline**: Converts a FASTA file with multiline sequences into a file where each sequence is written on a single line.

   **Arguments**:
   - **input_fasta** — path to the input FASTA file.
   - **output_fasta** — path to the output FASTA file (optional).

   **Example usage**:
   ```bash
   python bio_files_processor.py convert_fasta input.fasta --output_fasta output.fasta
   ```

2. **parse_blast_output**: Parses a BLAST output file and extracts the descriptions of the best matches.

   **Arguments**:
   - **input_file** — path to the BLAST output file (in txt format).
   - **output_file** — path to the output file where descriptions will be saved.

   **Example usage**:
   ```bash
   python bio_files_processor.py parse_blast input_blast.txt output_descriptions.txt
   ```

3. **select_genes_from_gbk_to_fasta**: Extracts protein sequences from a GBK file for specified genes and their neighbors, and saves them to a FASTA file.

   **Arguments**:
   - **input_gbk** — path to the input GBK file.
   - **genes** — list of target genes.
   - **n_before** — number of genes to include before the target gene.
   - **n_after** — number of genes to include after the target gene.
   - **output_fasta** — name of the output FASTA file.

   **Example usage**:
   ```bash
   python bio_files_processor.py select_genes input.gbk geneA geneB --n_before 1 --n_after 1 --output_fasta output.fasta
   ```

### 4. filter_fastq (legacy)

If you want to use the function with a dictionary of sequences rather than a FASTQ file, you can refer to the older `filter_fastq` function. This function works with a dictionary format for FASTQ sequences (name, (sequence, quality)).

#### Arguments:
- **seqs** — a dictionary with FASTQ sequences (key: name, value: tuple of sequence and quality).
- **gc_bounds** — GC content bounds for filtering (default is (0, 100)).
- **length_bounds** — sequence length bounds for filtering (default is (0, 2**32)).
- **quality_threshold** — threshold for the average quality score (default is 0).

#### Example usage:
```python
from main_script import filter_fastq

seqs = {
    "seq1": ("ATGC", "IIII"),
    "seq2": ("CGTGA", "HHHHH")
}
gc_bounds = (40, 60)
length_bounds = (4, 10)
quality_threshold = 30

filtered_seqs = filter_fastq(seqs, gc_bounds, length_bounds, quality_threshold)
print(filtered_seqs)
```

## Project Structure

```
bioinf-utils/
├─ README.md
├─ main_script.py
├─ utils/
│  ├── __init__.py
│  ├── fastq_module.py
│  ├── dna_rna_tools_test.py
│  ├── fastq_filtrator_test.py
│  └── example_data.py
├─ bio_files_processor.py
```

## Screenshots
Here is a screenshot of test results:
```