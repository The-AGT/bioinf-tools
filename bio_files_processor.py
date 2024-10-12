import argparse
import os
import re
from typing import List, Optional, Dict


def convert_multiline_fasta_to_oneline(
    input_fasta: str,
    output_fasta: Optional[str] = None
):
    """
    Converts a FASTA file with multiline sequences into a file where
    each sequence is written on a single line.

    :param input_fasta: Path to the input FASTA file.
    :param output_fasta: Path to the output FASTA file. If not specified,
                         the output file is created automatically.
    """
    if output_fasta is None:
        output_fasta = f"{os.path.splitext(input_fasta)[0]}_oneline.fasta"

    with open(input_fasta, 'r') as infile, open(output_fasta, 'w') as outfile:
        sequence_lines = []
        header = None

        for line in infile:
            line = line.rstrip()
            if line.startswith('>'):
                if header and sequence_lines:
                    outfile.write(f"{header}\n{''.join(sequence_lines)}\n")
                header = line
                sequence_lines = []
            else:
                sequence_lines.append(line)

        if header and sequence_lines:
            outfile.write(f"{header}\n{''.join(sequence_lines)}\n")


def parse_blast_output(
    input_file: str,
    output_file: Optional[str] = None
):
    """
    Parses BLAST output data and extracts the best match descriptions.

    :param input_file: Path to the BLAST output file (txt format).
    :param output_file: Path to the output file where descriptions
                        will be saved. If not specified, the output file
                        is created automatically.
    """
    if output_file is None:
        output_file = f"{os.path.splitext(input_file)[0]}_parsed.txt"

    with open(input_file, 'r') as infile, open("temp_file.txt", 'w') as tempfile:
        parsing = False

        for line in infile:
            if "Sequences producing significant alignments:" in line:
                parsing = True
                continue

            if parsing:
                next_line = infile.readline().strip()
                next_line = infile.readline().strip()

                if next_line:
                    description = next_line.split('...')[0].strip()
                    description = re.split(r"\s{2,}", description)[0]
                    tempfile.write(description + "\n")
                    parsing = False

    sort_and_remove_duplicates(
        temp_file="temp_file.txt",
        output_file=output_file
    )


def sort_and_remove_duplicates(
    temp_file: str,
    output_file: str
):
    """
    Sorts and removes duplicate lines from the temporary file and writes
    the result to the output file.

    :param temp_file: Path to the temporary file used for storing
                      unsorted descriptions.
    :param output_file: Path to the output file where the sorted, unique
                        descriptions will be saved.
    """
    lines = set()

    with open(temp_file, 'r') as tempfile:
        for line in tempfile:
            lines.add(line.strip())

    sorted_lines = sorted(lines)

    with open(output_file, 'w') as outfile:
        for line in sorted_lines:
            outfile.write(line + "\n")

    os.remove(temp_file)


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
    :param n_before: Number of genes before the target gene.
                     Default is 1.
    :param n_after: Number of genes after the target gene.
                    Default is 1.
    :param output_fasta: Name of the output FASTA file.
                         Default is 'output.fasta'.
    """
    def parse_features(lines: List[str]) -> List[Dict]:
        """Parses the FEATURES section in a GBK file."""
        features = []
        current_feature = None
        in_feature_section = False
        current_key = None
        current_value = ''

        feature_start_re = re.compile(r'^\s{5}\S')
        qualifier_start_re = re.compile(r'^\s{21}\S')

        for line in lines:
            if line.startswith('FEATURES'):
                in_feature_section = True
                continue
            if not in_feature_section:
                continue
            if line.startswith('ORIGIN'):
                if current_feature:
                    if current_key:
                        current_feature['qualifiers'][current_key] = current_value.strip('"')
                    features.append(current_feature)
                break
            if feature_start_re.match(line):
                if current_feature:
                    if current_key:
                        current_feature['qualifiers'][current_key] = current_value.strip('"')
                    features.append(current_feature)
                parts = line.strip().split()
                if len(parts) >= 2:
                    feature_type = parts[0]
                    location = ' '.join(parts[1:])
                    current_feature = {'type': feature_type,
                                       'location': location,
                                       'qualifiers': {}}
                else:
                    current_feature = None
                current_key = None
                current_value = ''
            elif qualifier_start_re.match(line):
                qualifier_line = line.strip()
                if qualifier_line.startswith('/'):
                    if current_key:
                        current_feature['qualifiers'][current_key] = current_value.strip('"')
                    if '=' in qualifier_line:
                        key, value = qualifier_line.lstrip('/').split('=', 1)
                        current_key = key.strip()
                        current_value = value.strip('"')
                        if current_value.endswith('"') and current_value.count('"') % 2 == 1:
                            current_value = current_value[:-1]
                    else:
                        current_key = qualifier_line.lstrip('/')
                        current_value = ''
                else:
                    current_value += qualifier_line.strip()
        if current_feature:
            if current_key:
                current_feature['qualifiers'][current_key] = current_value.strip('"')
            features.append(current_feature)

        return features

    def find_target_genes(
        features: List[Dict],
        target_genes: List[str]
    ) -> List[Dict]:
        """Finds the target genes and surrounding genes."""
        gene_indices = []
        for idx, feature in enumerate(features):
            if feature['type'] == 'CDS':
                gene_names = [
                    feature['qualifiers'].get(qualifier, '')
                    for qualifier in ['gene', 'locus_tag', 'product']
                ]
                if any(any(target_gene in gene_name for gene_name in gene_names)
                       for target_gene in target_genes):
                    gene_indices.append(idx)

        selected_genes = []
        for idx in gene_indices:
            start_idx = max(0, idx - n_before)
            end_idx = min(len(features), idx + n_after + 1)
            selected_genes.extend(features[start_idx:end_idx])

        seen = set()
        unique_genes = []
        for gene in selected_genes:
            gene_id = id(gene)
            if gene_id not in seen:
                unique_genes.append(gene)
                seen.add(gene_id)

        return unique_genes

    def write_fasta(
        genes: List[Dict],
        output_fasta: str
    ):
        """Writes the selected genes to the output FASTA file."""
        with open(output_fasta, 'w') as outfile:
            for gene in genes:
                qualifiers = gene['qualifiers']
                protein_seq = qualifiers.get('translation', '').replace('\n', '').replace(' ', '')
                if protein_seq:
                    gene_name = (qualifiers.get('gene')
                                 or qualifiers.get('locus_tag')
                                 or qualifiers.get('product')
                                 or 'unknown_gene')
                    outfile.write(f">{gene_name}\n{protein_seq}\n")

    if not os.path.exists(input_gbk):
        raise FileNotFoundError(f"File {input_gbk} not found.")

    with open(input_gbk, 'r') as file:
        lines = file.readlines()

    if not any("FEATURES" in line for line in lines):
        raise ValueError("The GBK file does not contain a FEATURES section.")

    features = parse_features(lines)
    selected_genes = find_target_genes(features, genes)

    if not selected_genes:
        raise ValueError("No target genes were found in the GBK file.")

    write_fasta(selected_genes, output_fasta)


def main():
    """Main function to parse command-line arguments and call appropriate functions."""
    parser = argparse.ArgumentParser(
        description="Utility for processing bioinformatics files using standard Python libraries."
    )
    subparsers = parser.add_subparsers(dest="command", help="Available commands")

    parser_fasta = subparsers.add_parser(
        'convert_fasta',
        help='Convert a FASTA file with multiline sequences to single-line sequences.'
    )
    parser_fasta.add_argument('input_fasta', help='Path to the input FASTA file.')
    parser_fasta.add_argument(
        '--output_fasta',
        help='Path to the output FASTA file. If not specified, it will be created automatically.',
        default=None
    )

    parser_blast = subparsers.add_parser(
        'parse_blast',
        help='Parse BLAST output data and extract the best match descriptions.'
    )
    parser_blast.add_argument('input_file', help='Path to the BLAST output file (txt format).')
    parser_blast.add_argument('output_file', help='Path to the output file where descriptions will be saved.')

    parser_gbk = subparsers.add_parser(
        'select_genes',
        help='Extract genes from a GBK file and save them to a FASTA file.'
    )
    parser_gbk.add_argument('input_gbk', help='Path to the input GBK file.')
    parser_gbk.add_argument('genes', nargs='+', help='List of target genes separated by space.')
    parser_gbk.add_argument('--n_before', type=int, default=1, help='Number of genes before the target gene. Default is 1.')
    parser_gbk.add_argument('--n_after', type=int, default=1, help='Number of genes after the target gene. Default is 1.')
    parser_gbk.add_argument('--output_fasta', default='output.fasta', help='Name of the output FASTA file. Default is output.fasta.')

    args = parser.parse_args()

    if args.command == 'convert_fasta':
        convert_multiline_fasta_to_oneline(args.input_fasta, args.output_fasta)
    elif args.command == 'parse_blast':
        parse_blast_output(args.input_file, args.output_file)
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
