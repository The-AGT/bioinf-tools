#!/usr/bin/env python3
"""
Module for biological sequence operations and FASTQ filtering.
Uses Biopython for FASTQ parsing and filtering.
"""

from abc import ABC, abstractmethod
from collections import Counter
import argparse
import os
from typing import Tuple, List, Union, Dict, Any

from Bio import SeqIO


class BiologicalSequence(ABC):
    """
    Abstract base class for biological sequences.

    Implements:
      - Length, indexing, and slicing (returning same type)
      - String representation (__str__, __repr__)
      - Equality and concatenation
      - Alphabet validation
    """

    def __init__(self, sequence: str):
        self._sequence = sequence
        if not self.check_alphabet():
            raise ValueError(
                f"Sequence contains invalid characters: {sequence}"
            )

    @property
    def sequence(self) -> str:
        """Return the underlying sequence string."""
        return self._sequence

    def __len__(self) -> int:
        return len(self._sequence)

    def __getitem__(self, index):
        result = self._sequence[index]
        if isinstance(index, slice):
            return type(self)(result)
        return result

    def __str__(self) -> str:
        return self._sequence

    def __repr__(self) -> str:
        return f"{self.__class__.__name__}({self._sequence})"

    def __eq__(self, other) -> bool:
        if isinstance(other, BiologicalSequence):
            return self._sequence == other._sequence
        return False

    def __add__(self, other):
        if isinstance(other, BiologicalSequence):
            return type(self)(self._sequence + other._sequence)
        if isinstance(other, str):
            return type(self)(self._sequence + other)
        raise TypeError("Can only concatenate BiologicalSequence or str")

    def check_alphabet(self) -> bool:
        """
        Check if all characters in the sequence are in the valid alphabet.
        """
        return all(base in self.valid_alphabet() for base in self._sequence)

    @abstractmethod
    def valid_alphabet(self) -> set:
        """
        Return a set of valid characters for the sequence.
        Must be implemented in subclasses.
        """
        pass


class NucleicAcidSequence(BiologicalSequence):
    """
    Abstract class for nucleic acid sequences (DNA and RNA).

    Implements reverse, complement, and reverse complement.
    """

    @abstractmethod
    def _complement_mapping(self) -> dict:
        """
        Return a mapping dictionary for the complement of each base.
        Must be implemented in subclasses.
        """
        pass

    def reverse(self):
        """
        Return a new object of the same type with the reversed sequence.
        """
        return type(self)(self._sequence[::-1])

    def complement(self):
        """
        Return a new object with the complementary sequence.
        """
        mapping = self._complement_mapping()
        try:
            comp_seq = "".join(mapping[base] for base in self._sequence)
        except KeyError as exc:
            raise ValueError(f"Invalid base encountered: {exc}")
        return type(self)(comp_seq)

    def reverse_complement(self):
        """
        Return a new object with the reverse complement sequence.
        """
        return self.complement().reverse()


class DNASequence(NucleicAcidSequence):
    """
    Class for DNA sequences.

    Implements:
      - Valid alphabet (ATCG, case-sensitive)
      - Complement mapping for DNA
      - Transcription to RNA (returns an RNASequence)
    """

    def valid_alphabet(self) -> set:
        return set("ATCGatcg")

    def _complement_mapping(self) -> dict:
        return {
            "A": "T",
            "T": "A",
            "C": "G",
            "G": "C",
            "a": "t",
            "t": "a",
            "c": "g",
            "g": "c",
        }

    def transcribe(self):
        """
        Transcribe DNA to RNA by replacing T/t with U/u.
        Returns an RNASequence object.
        """
        rna_seq = self._sequence.replace("T", "U").replace("t", "u")
        return RNASequence(rna_seq)


class RNASequence(NucleicAcidSequence):
    """
    Class for RNA sequences.

    Implements:
      - Valid alphabet (AUCG, case-sensitive)
      - Complement mapping for RNA
    """

    def valid_alphabet(self) -> set:
        return set("AUCGaucg")

    def _complement_mapping(self) -> dict:
        return {
            "A": "U",
            "U": "A",
            "C": "G",
            "G": "C",
            "a": "u",
            "u": "a",
            "c": "g",
            "g": "c",
        }


class AminoAcidSequence(BiologicalSequence):
    """
    Class for amino acid sequences.

    Implements:
      - Valid alphabet (20 standard amino acids, case-sensitive)
      - Amino acid composition as a dictionary
      - Hydrophobicity score using the Kyte-Doolittle scale
    """

    def valid_alphabet(self) -> set:
        return set("ACDEFGHIKLMNPQRSTVWY"
                   "acdefghiklmnpqrstvwy")

    def amino_acid_composition(self) -> dict:
        """
        Return a dictionary with the count of each amino acid.
        """
        return dict(Counter(self._sequence))

    def hydrophobicity_score(self) -> float:
        """
        Compute the total hydrophobicity using the Kyte-Doolittle scale.
        """
        kd = {
            "A": 1.8,
            "C": 2.5,
            "D": -3.5,
            "E": -3.5,
            "F": 2.8,
            "G": -0.4,
            "H": -3.2,
            "I": 4.5,
            "K": -3.9,
            "L": 3.8,
            "M": 1.9,
            "N": -3.5,
            "P": -1.6,
            "Q": -3.5,
            "R": -4.5,
            "S": -0.8,
            "T": -0.7,
            "V": 4.2,
            "W": -0.9,
            "Y": -1.3,
        }
        return sum(kd.get(aa.upper(), 0) for aa in self._sequence)


def filter_fastq(input_fastq: Any,
                 output_fastq: Union[str, None] = None,
                 gc_bounds: Union[Tuple[float, float], float, int] =
                 (0.0, 1.0),
                 length_bounds: Union[Tuple[int, int], int] =
                 (0, 2 ** 32),
                 quality_threshold: int = 0
                 ) -> Union[Dict[str, Any], None]:
    """
    Filter a FASTQ file or in-memory FASTQ records using Biopython.

    Filtering criteria:
      - Sequence length within length_bounds.
      - Average quality (computed from Phred scores) at least
        quality_threshold.
      - GC fraction (computed from the sequence) within gc_bounds.
        If gc_bounds values are greater than 1, they are assumed to be
        percentages and are converted to fractions.

    When input_fastq is a file path (str), the file is parsed via SeqIO.
    When input_fastq is a dict, it is assumed to map record IDs to tuples
    of the form (sequence, quality, extra_info). In this case the original
    tuple is returned if the record passes filtering.

    If output_fastq is provided (as a file path), the filtered FASTQ
    records are written to that file and None is returned.
    Otherwise, a dict of filtered records is returned.

    :param input_fastq: Path to the FASTQ file or a dict of FASTQ records.
    :param output_fastq: Optional output file path.
    :param gc_bounds: GC fraction bounds (or a single number for strict
                      equality). If >1, assumed percentage.
    :param length_bounds: Tuple (or int) for sequence length bounds.
    :param quality_threshold: Minimum average quality.
    :raises FileExistsError: If output file exists.
    :return: A dict of filtered records or None if written to file.
    """
    # If gc_bounds is a single number, convert it to a tuple.
    if isinstance(gc_bounds, (int, float)):
        gc_bounds = (float(gc_bounds), float(gc_bounds))

    # If gc_bounds are given as percentages, convert them.
    if (isinstance(gc_bounds, (list, tuple)) and len(gc_bounds) == 2 and
            gc_bounds[1] > 1):
        gc_bounds = (gc_bounds[0] / 100.0, gc_bounds[1] / 100.0)

    # If length_bounds is an int, convert it to a tuple.
    if isinstance(length_bounds, int):
        length_bounds = (0, length_bounds)

    # If input_fastq is a file path (string).
    if isinstance(input_fastq, str):
        filtered_records = []
        for record in SeqIO.parse(input_fastq, "fastq"):
            seq = str(record.seq)
            seq_len = len(seq)
            if not (length_bounds[0] <= seq_len <= length_bounds[1]):
                continue
            qualities = record.letter_annotations.get("phred_quality", [])
            if not qualities:
                continue
            avg_quality = sum(qualities) / seq_len
            if avg_quality < quality_threshold:
                continue
            # Compute gc_fraction manually.
            gc_fraction = ((seq.upper().count("G") +
                            seq.upper().count("C")) / seq_len)
            if not (gc_bounds[0] <= gc_fraction <= gc_bounds[1]):
                continue
            filtered_records.append(record)
        if output_fastq:
            if os.path.exists(output_fastq):
                raise FileExistsError(
                    f"Output file '{output_fastq}' already exists."
                )
            with open(output_fastq, "w") as out_handle:
                SeqIO.write(filtered_records, out_handle, "fastq")
            return None
        else:
            result = {}
            for record in filtered_records:
                key = (record.id if record.id.startswith("@")
                       else "@" + record.id)
                qual_list = record.letter_annotations.get("phred_quality", [])
                quality_str = "".join(chr(q + 33) for q in qual_list)
                result[key] = (str(record.seq), quality_str, quality_str)
            return result
    # If input_fastq is a dict of FASTQ records.
    elif isinstance(input_fastq, dict):
        result = {}
        for key, value in input_fastq.items():
            if (not isinstance(value, (list, tuple)) or
                    len(value) < 2):
                continue
            seq = value[0]
            quality = value[1]
            seq_len = len(seq)
            if seq_len == 0 or not (length_bounds[0] <= seq_len <=
                                     length_bounds[1]):
                continue
            avg_quality = (sum(ord(c) - 33 for c in quality) /
                           len(quality))
            if avg_quality < quality_threshold:
                continue
            gc_fraction = ((seq.upper().count("G") +
                            seq.upper().count("C")) / seq_len)
            if not (gc_bounds[0] <= gc_fraction <= gc_bounds[1]):
                continue
            result[key] = value
        return result
    else:
        raise TypeError("input_fastq must be a file path or a dict of records")


def run_sequence_operation(seq_type: str, sequences: List[str],
                           operation: str) -> Union[str, List]:
    """
    Run the specified operation on given sequences.

    Supported operations:
      DNA: reverse, complement, reverse_complement, transcribe
      RNA: reverse, complement, reverse_complement
      protein: reverse, amino_acid_composition,
               hydrophobicity_score

    If multiple sequences are provided, returns a list.
    """
    result = []
    for seq in sequences:
        if seq_type == "dna":
            obj = DNASequence(seq)
            if operation == "reverse":
                res = obj.reverse()
            elif operation == "complement":
                res = obj.complement()
            elif operation == "reverse_complement":
                res = obj.reverse_complement()
            elif operation == "transcribe":
                res = obj.transcribe()
            else:
                raise ValueError(
                    f"Invalid operation '{operation}' for DNA sequence."
                )
        elif seq_type == "rna":
            obj = RNASequence(seq)
            if operation == "reverse":
                res = obj.reverse()
            elif operation == "complement":
                res = obj.complement()
            elif operation == "reverse_complement":
                res = obj.reverse_complement()
            else:
                raise ValueError(
                    f"Invalid operation '{operation}' for RNA sequence."
                )
        elif seq_type == "protein":
            obj = AminoAcidSequence(seq)
            if operation == "reverse":
                res = obj.reverse()
            elif operation == "amino_acid_composition":
                res = obj.amino_acid_composition()
            elif operation == "hydrophobicity_score":
                res = obj.hydrophobicity_score()
            else:
                raise ValueError(
                    f"Invalid operation '{operation}' for protein sequence."
                )
        else:
            raise ValueError(f"Unknown sequence type: {seq_type}")
        result.append(res)
    return result if len(result) > 1 else result[0]


def main():
    """
    Main function for the command-line interface.

    Supports two subcommands:
      - filter: FASTQ file filtering using Biopython.
      - run: Operations on sequences (DNA, RNA, protein).
    """
    parser = argparse.ArgumentParser(
        description=("Operations on biological sequences and FASTQ filtering "
                     "using Biopython.")
    )
    subparsers = parser.add_subparsers(
        dest="command", help="Available commands"
    )

    # Subcommand for FASTQ filtering.
    parser_filter = subparsers.add_parser(
        "filter", help="Filter a FASTQ file"
    )
    parser_filter.add_argument(
        "input_fastq",
        help="Path to the input FASTQ file or dict of records"
    )
    parser_filter.add_argument(
        "output_fastq", nargs="?",
        help=("Path to the output FASTQ file (if not provided, returns dict)")
    )
    parser_filter.add_argument(
        "--gc_bounds", type=float, nargs="*",
        default=(0.0, 1.0),
        help=("GC fraction bounds: one or two numbers "
              "(if >1, assumed percentage)")
    )
    parser_filter.add_argument(
        "--length_bounds", type=int, nargs="*",
        default=(0, 2 ** 32),
        help=("Sequence length bounds: one or two numbers "
              "(lower and upper bounds)")
    )
    parser_filter.add_argument(
        "--quality_threshold", type=int, default=0,
        help="Minimum average quality"
    )

    # Subcommand for sequence operations.
    parser_run = subparsers.add_parser(
        "run", help="Run operations on sequences"
    )
    parser_run.add_argument(
        "seq_type", choices=["dna", "rna", "protein"],
        help="Sequence type: dna, rna or protein"
    )
    parser_run.add_argument(
        "sequences", nargs="+",
        help=("Sequences followed by the operation "
              "(last argument is the operation)")
    )

    args = parser.parse_args()

    if args.command == "filter":
        gc_bounds = tuple(args.gc_bounds)
        if len(gc_bounds) == 1:
            gc_bounds = (0.0, gc_bounds[0])
        elif len(gc_bounds) != 2:
            gc_bounds = (0.0, 1.0)

        length_bounds = tuple(args.length_bounds)
        if len(length_bounds) == 1:
            length_bounds = (0, length_bounds[0])
        elif len(length_bounds) != 2:
            length_bounds = (0, 2 ** 32)

        try:
            result = filter_fastq(
                args.input_fastq,
                args.output_fastq,
                gc_bounds=gc_bounds,
                length_bounds=length_bounds,
                quality_threshold=args.quality_threshold
            )
            if result is not None:
                print(result)
        except FileExistsError as exc:
            print(exc)
            exit(1)
    elif args.command == "run":
        *seqs, op = args.sequences
        try:
            result = run_sequence_operation(args.seq_type, seqs, op)
            print(result)
        except ValueError as exc:
            print(exc)
            exit(1)
    else:
        parser.print_help()


if __name__ == "__main__":
    main()
