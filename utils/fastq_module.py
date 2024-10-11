# fastq_module.py

from typing import Tuple, Generator
import os


def read_fastq(input_fastq: str) -> Generator[Tuple[str, Tuple[str, str]],
                                              None, None]:
    """
    Generator for reading a FASTQ file and yielding sequences one by one
    as a dictionary-like structure.

    :param input_fastq: Path to the input FASTQ file.
    :yield: A tuple (name, (sequence, quality)) representing the FASTQ record.
    """
    with open(input_fastq, 'r') as file:
        while True:
            header = file.readline().rstrip()
            if not header:
                break  # End of file
            sequence = file.readline().rstrip()
            file.readline().rstrip()  # Skip separator line ('+')
            quality = file.readline().rstrip()
            if not quality:
                break  # Incomplete record
            name = header[1:]  # Remove leading '@'
            yield name, (sequence, quality)


def write_fastq(name: str, sequence: str, quality: str, output_handle):
    """
    Writes a single FASTQ sequence to the output file.

    :param name: Name of the sequence.
    :param sequence: Nucleotide sequence.
    :param quality: Quality string.
    :param output_handle: File handle for writing.
    """
    output_handle.write(f"@{name}\n{sequence}\n+\n{quality}\n")


def is_within_bounds(value: float, bounds: Tuple[int, int]) -> bool:
    """
    Checks if a value is within the specified bounds.

    :param value: The value to check.
    :param bounds: A tuple containing the lower and upper bounds.
    :return: True if the value is within bounds, otherwise False.
    """
    lower, upper = bounds
    return lower <= value <= upper


def is_sequence_valid(
    sequence: str,
    quality: str,
    gc_bounds: Tuple[int, int] = (0, 100),
    length_bounds: Tuple[int, int] = (0, 2**32),
    quality_threshold: int = 0
) -> bool:
    """
    Validates if a FASTQ sequence meets the filtering criteria.

    :param sequence: Nucleotide sequence.
    :param quality: Quality string.
    :param gc_bounds: GC content bounds (%).
    :param length_bounds: Sequence length bounds.
    :param quality_threshold: Minimum average quality score.
    :return: True if the sequence passes the filters, otherwise False.
    """
    # Check sequence length
    seq_length = len(sequence)
    if not is_within_bounds(seq_length, length_bounds):
        return False

    # Calculate GC content
    gc_count = sequence.upper().count('G') + sequence.upper().count('C')
    gc_content_value = (gc_count / seq_length) * 100 if seq_length > 0 else 0.0
    if not is_within_bounds(gc_content_value, gc_bounds):
        return False

    # Calculate average quality
    total_quality = sum(ord(char) - 33 for char in quality)
    avg_quality = total_quality / len(quality) if len(quality) > 0 else 0.0
    if avg_quality < quality_threshold:
        return False

    return True


def filter_fastq_file(
    input_fastq: str,  # Path to the input FASTQ file
    output_fastq: str,  # Path to the output FASTQ file
    gc_bounds: Tuple[int, int] = (0, 100),
    length_bounds: Tuple[int, int] = (0, 2**32),
    quality_threshold: int = 0
):
    """
    Reads sequences from an input FASTQ file,
    filters them based on GC content, length, and quality,
    and writes the filtered sequences to an output FASTQ file.

    :param input_fastq: Path to the input FASTQ file.
    :param output_fastq: Name of the output FASTQ file.
    :param gc_bounds: GC content bounds (%).
    :param length_bounds: Sequence length bounds.
    :param quality_threshold: Minimum average quality score.
    :raises FileExistsError: If the output file already exists.
    """
    output_dir = 'filtered'
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    output_path = os.path.join(output_dir, output_fastq)

    if os.path.exists(output_path):
        raise FileExistsError(f"Output file '{output_path}' already exists.")

    # Process sequences one by one
    with open(output_path, 'w') as output_handle:
        for name, (sequence, quality) in read_fastq(input_fastq):
            if is_sequence_valid(sequence, quality, gc_bounds,
                                 length_bounds, quality_threshold):
                write_fastq(name, sequence, quality, output_handle)
