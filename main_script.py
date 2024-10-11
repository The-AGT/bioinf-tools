# main_script.py

from typing import List, Union
from utils.helpers import (
    is_valid_sequence, transcribe, reverse,
    complement, reverse_complement, gc_content
)
from utils.fastq_module import filter_fastq_file


def run_dna_rna_tools(*args: str) -> Union[str, List[str]]:
    """
    Runs specified operations on DNA or RNA sequences.
    :param args: Sequences (DNA or RNA) and the operation.
                 The operation must be the last argument.
    :return: The result of applying the operation to the sequence(s).
    :raises ValueError: If no sequences or operation are provided.
    """
    if len(args) < 2:
        raise ValueError(
            "At least one sequence and one operation must be provided."
        )

    # The last argument is the operation
    *sequences, operation = args

    valid_bases = {'A', 'T', 'C', 'G', 'U', 'a', 't', 'c', 'g', 'u'}
    valid_procedures = {'transcribe', 'reverse', 'complement',
                        'reverse_complement', 'gc_content'}

    if operation not in valid_procedures:
        raise ValueError(f"Invalid operation: {operation}")

    results = []

    for seq in sequences:
        if not is_valid_sequence(seq, valid_bases):
            raise ValueError(f"Invalid sequence: {seq}")
        if operation == 'transcribe':
            results.append(transcribe(seq))
        elif operation == 'reverse':
            results.append(reverse(seq))
        elif operation == 'complement':
            results.append(complement(seq))
        elif operation == 'reverse_complement':
            results.append(reverse_complement(seq))
        elif operation == 'gc_content':
            results.append(gc_content(seq))

    return results if len(results) > 1 else results[0]


# Alias for run_dna_rna_tools to match tests
rdrt = run_dna_rna_tools


def main():
    """
    Main function for processing command-line arguments
    and running FASTQ filtering or DNA/RNA operations.
    """
    import argparse

    parser = argparse.ArgumentParser(
        description='Filter FASTQ files or run DNA/RNA sequence operations.'
    )

    subparsers = parser.add_subparsers(dest="command",
                                       help="Available commands")

    # Subcommand for filtering FASTQ files
    parser_filter = subparsers.add_parser(
        "filter", help="Filter a FASTQ file"
    )
    parser_filter.add_argument(
        'input_fastq',
        help='Path to the input FASTQ file.'
    )
    parser_filter.add_argument(
        'output_fastq',
        help='Name of the output FASTQ file.'
    )
    parser_filter.add_argument(
        '--gc_bounds',
        type=float,
        nargs='*',
        default=(0, 100),
        help='GC content bounds: one or two numbers (lower and upper).'
    )
    parser_filter.add_argument(
        '--length_bounds',
        type=int,
        nargs='*',
        default=(0, 2**32),
        help='Sequence length bounds: one or two numbers (lower and upper).'
    )
    parser_filter.add_argument(
        '--quality_threshold',
        type=int,
        default=0,
        help='Minimum average quality score.'
    )

    # Subcommand for running DNA/RNA sequence operations
    parser_dna_rna = subparsers.add_parser(
        "run", help="Run operations on DNA/RNA sequences"
    )
    parser_dna_rna.add_argument(
        'sequences',
        nargs='+',
        help="Sequences and operation (the last argument is the operation)"
    )

    args = parser.parse_args()

    if args.command == "filter":
        # Process gc_bounds and length_bounds for FASTQ filtering
        gc_bounds = tuple(args.gc_bounds)
        if len(gc_bounds) == 1:
            gc_bounds = (0, gc_bounds[0])
        elif len(gc_bounds) == 2:
            pass
        else:
            gc_bounds = (0, 100)

        length_bounds = tuple(args.length_bounds)
        if len(length_bounds) == 1:
            length_bounds = (0, length_bounds[0])
        elif len(length_bounds) == 2:
            pass
        else:
            length_bounds = (0, 2**32)

        try:
            filter_fastq_file(
                args.input_fastq,
                args.output_fastq,
                gc_bounds=gc_bounds,
                length_bounds=length_bounds,
                quality_threshold=args.quality_threshold
            )
        except FileExistsError as e:
            print(e)
            exit(1)

    elif args.command == "run":
        # Execute operation on sequences (run_dna_rna_tools)
        try:
            result = run_dna_rna_tools(*args.sequences)
            print(result)
        except ValueError as e:
            print(e)
            exit(1)

    else:
        parser.print_help()


if __name__ == '__main__':
    main()
