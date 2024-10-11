# bio_files_processor.py

import argparse
import os
from typing import List, Dict, Optional


def convert_multiline_fasta_to_oneline(input_fasta: str,
                                       output_fasta: Optional[str] = None):
    """
    Converts a FASTA file with multiline sequences into a file where
    each sequence is written on a single line.

    :param input_fasta: Path to the input FASTA file.
    :param output_fasta: Path to the output FASTA file. If not specified,
                         an output file is created automatically.
    """
    if output_fasta is None:
        output_fasta = f"{os.path.splitext(input_fasta)[0]}_oneline.fasta"

    with open(input_fasta, 'r') as infile, open(output_fasta, 'w') as outfile:
        sequence_lines = []
        header = None

        for line in infile:
            line = line.rstrip()
            if line.startswith('>'):
                if header:
                    # Write the current sequence to the file
                    outfile.write(f"{header}\n{''.join(sequence_lines)}\n")
                header = line  # Store new header
                sequence_lines = []  # Reset sequence
            else:
                sequence_lines.append(line)

        # Write the last sequence
        if header:
            outfile.write(f"{header}\n{''.join(sequence_lines)}\n")


def parse_blast_output(input_file: str, output_file: str):
    """
    Parses BLAST output data and extracts the best match descriptions.

    :param input_file: Path to the BLAST output file (txt format).
    :param output_file: Path to the output file where descriptions
                        will be saved.
    """
    descriptions = set()

    with open(input_file, 'r') as infile:
        lines = infile.readlines()
        inside_query = False
        for idx, line in enumerate(lines):
            line = line.strip()
            if line.startswith('Query='):
                inside_query = True
            elif inside_query and 'Sequences producing significant alignments:' in line:  # noqa: E501
                # Extract the description of the first match after this header
                description_line = lines[idx + 2].strip()
                if description_line:
                    description = ' '.join(description_line.split()[1:])
                    descriptions.add(description)
                inside_query = False

    with open(output_file, 'w') as outfile:
        for desc in sorted(descriptions):
            outfile.write(f"{desc}\n")


def select_genes_from_gbk_to_fasta(
    input_gbk: str,
    genes: List[str],
    n_before: int = 1,
    n_after: int = 1,
    output_fasta: str = 'output.fasta'
):
    """
    Extracts protein sequences from a GBK file for specified genes
    and their neighbors, and saves them to a FASTA file.

    :param input_gbk: Path to the input GBK file.
    :param genes: List of target genes.
    :param n_before: Number of genes before the target gene. Default is 1.
    :param n_after: Number of genes after the target gene. Default is 1.
    :param output_fasta: Name of the output FASTA file.
                         Default is 'output.fasta'.
    """
    def parse_features(lines: List[str]) -> List[Dict]:
        """Parses the FEATURES section of a GBK file."""
        features = []
        current_feature = None
        in_feature_section = False

        for line in lines:
            if line.startswith('FEATURES'):
                in_feature_section = True
            elif in_feature_section and line.startswith('     '):
                # Start of a new feature
                if current_feature:
                    features.append(current_feature)
                current_feature = {
                    'type': line.split()[0],
                    'location': line.split()[1:],
                    'qualifiers': {}
                }
            elif current_feature and line.startswith('                     '):
                # Process qualifiers
                if '=' in line:
                    key, value = line.strip().lstrip('/').split('=', 1)
                    current_feature['qualifiers'][key.strip()] = value.strip('"')  # noqa: E501

        if current_feature:
            features.append(current_feature)

        return features

    def find_target_genes(features: List[Dict],
                          target_genes: List[str]) -> List[Dict]:
        """Finds the target genes and surrounding genes."""
        gene_indices = [
            idx for idx, feature in enumerate(features)
            if (
                feature['type'] == 'CDS'
                and 'gene' in feature['qualifiers']
                and any(gene in feature['qualifiers']['gene']
                        for gene in target_genes)
            )
        ]

        selected_genes = []
        for idx in gene_indices:
            start_idx = max(0, idx - n_before)
            end_idx = min(len(features), idx + n_after + 1)
            selected_genes.extend(features[start_idx:end_idx])

        return selected_genes

    def write_fasta(genes: List[Dict], output_fasta: str):
        """Writes the selected genes to the output FASTA file."""
        with open(output_fasta, 'w') as outfile:
            for gene in genes:
                qualifiers = gene['qualifiers']
                protein_seq = qualifiers.get('translation', '').replace('\n', '').replace(' ', '')  # noqa: E501
                if protein_seq:
                    gene_name = qualifiers.get('gene', 'unknown_gene')
                    outfile.write(f">{gene_name}\n{protein_seq}\n")

    # Read the GBK file
    with open(input_gbk, 'r') as file:
        lines = file.readlines()

    # Parse the features
    features = parse_features(lines)

    # Find the target genes
    selected_genes = find_target_genes(features, genes)

    # Write to the output FASTA file
    write_fasta(selected_genes, output_fasta)


def main():
    """Main function to parse command-line arguments and call appropriate functions."""  # noqa: E501
    parser = argparse.ArgumentParser(
        description='Utility for processing bioinformatics files using standard Python libraries.'  # noqa: E501
    )
    subparsers = parser.add_subparsers(dest='command',
                                       help='Available commands')

    # Subparser for the convert_multiline_fasta_to_oneline function
    parser_fasta = subparsers.add_parser(
        'convert_fasta',
        help='Convert a FASTA file with multiline sequences to single-line sequences.'  # noqa: E501
    )
    parser_fasta.add_argument(
        'input_fasta',
        help='Path to the input FASTA file.'
    )
    parser_fasta.add_argument(
        '--output_fasta',
        help='Path to the output FASTA file. If not specified, it will be created automatically.',  # noqa: E501
        default=None
    )

    # Subparser for the parse_blast_output function
    parser_blast = subparsers.add_parser(
        'parse_blast',
        help='Parse BLAST output data and extract the best match descriptions.'
    )
    parser_blast.add_argument(
        'input_file',
        help='Path to the BLAST output file (txt format).'
    )
    parser_blast.add_argument(
        'output_file',
        help='Path to the output file where descriptions will be saved.'
    )

    # Subparser for the select_genes_from_gbk_to_fasta function
    parser_gbk = subparsers.add_parser(
        'select_genes',
        help='Extract genes from a GBK file and save them to a FASTA file.'
    )
    parser_gbk.add_argument(
        'input_gbk',
        help='Path to the input GBK file.'
    )
    parser_gbk.add_argument(
        'genes',
        nargs='+',
        help='List of target genes separated by space.'
    )
    parser_gbk.add_argument(
        '--n_before',
        type=int,
        default=1,
        help='Number of genes before the target gene. Default is 1.'
    )
    parser_gbk.add_argument(
        '--n_after',
        type=int,
        default=1,
        help='Number of genes after the target gene. Default is 1.'
    )
    parser_gbk.add_argument(
        '--output_fasta',
        default='output.fasta',
        help='Name of the output FASTA file. Default is output.fasta.'
    )

    args = parser.parse_args()

    if args.command == 'convert_fasta':
        convert_multiline_fasta_to_oneline(
            args.input_fasta,
            args.output_fasta
        )
    elif args.command == 'parse_blast':
        parse_blast_output(
            args.input_file,
            args.output_file
        )
    elif args.command == 'select_genes':
        select_genes_from_gbk_to_fasta(
            args.input_gbk,
            args.genes,
            n_before=args.n_before,
            n_after=args.n_after,
            output_fasta=args.output_fasta
        )
    else:
        parser.print_help()


if __name__ == '__main__':
    main()
