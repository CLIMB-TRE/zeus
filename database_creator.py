import argparse
import os
import pybktree
import sys
import Levenshtein
import hashlib
import pyfastx
from glob import glob
import csv
import tqdm
import re


def md5_checksum(sequence: str) -> str:
    """Calculate the MD5 checksum of a given sequence.

    Args:
        sequence (str): The input sequence.

    Returns:
        str: The MD5 checksum of the sequence.
    """

    md5 = hashlib.md5()
    md5.update(sequence.encode("utf-8"))
    return md5.hexdigest()


def run(args: argparse.Namespace):

    names_dmp = os.path.join(args.taxonomy, "names.dmp")

    tax_id_name_lookup = {}
    name_tax_id_lookup = {}

    print(f"Loading NCBI taxonomy names from {names_dmp}", file=sys.stderr)
    with open(names_dmp, "r") as f:
        for line in f:
            fields = [i.lstrip() for i in line.split("\t|")]
            taxon_id, name, name_type = fields[0], fields[1], fields[3]
            if "scientific name" in name_type:
                tax_id_name_lookup[taxon_id] = name
                name_tax_id_lookup[name] = taxon_id

    print(
        f"Parsed {len(tax_id_name_lookup)} taxa from NCBI taxonomy names.dmp",
        file=sys.stderr,
    )

    name_tree = pybktree.BKTree(Levenshtein.distance, tax_id_name_lookup.values())

    print("Built BK-tree for taxon name lookup", file=sys.stderr)

    assemblies = glob(args.assembly_wildcard)

    print(f"Found {len(assemblies)} assemblies to process", file=sys.stderr)

    metadata_writer = csv.DictWriter(
        open(args.output_tsv, "wt", newline="\n", buffering=1),
        delimiter="\t",
        fieldnames=[
            "unique_accession",
            "accession_description",
            "taxon_id",
            "human_readable",
            "source",
            "source_accession",
            "sequence_md5",
            "sequence_length",
        ],
    )
    metadata_writer.writeheader()

    fasta_writer = open(args.output_fasta, "wt")

    if args.unmatched_assemblies:
        unmatched_assembly_writer = csv.DictWriter(
            open(args.unmatched_assemblies, "wt", newline="\n", buffering=1),
            delimiter="\t",
            fieldnames=["accession", "organism_name"],
        )
        unmatched_assembly_writer.writeheader()

    for assembly_dir in tqdm.tqdm(assemblies):

        metadata_record = {}

        accession = os.path.basename(assembly_dir)
        assembly_report_path = glob(
            os.path.join(assembly_dir, f"{accession}_*_assembly_report.txt")
        )[0]
        assembly_path = glob(os.path.join(assembly_dir, f"{accession}_*.fna.gz"))[0]

        with open(assembly_report_path, "r") as f:
            for line in f:
                if line.startswith("# Organism name:"):
                    organism_name = line.split(":")[1].strip()
                    break

        if organism_name not in tax_id_name_lookup:
            closest_matches = name_tree.find(re.sub(r"\(.*\)", "", organism_name), 5)
            if closest_matches:
                best_match = min(closest_matches, key=lambda x: x[0])
                matched_name = best_match[1]
                taxon_id = name_tax_id_lookup[matched_name]
                if args.verbose:
                    print(
                        f"Info: Organism name '{organism_name}' not found in taxonomy. "
                        f"Using closest match '{matched_name}' with taxon ID {taxon_id}.",
                        file=sys.stderr,
                    )
            else:
                if args.unmatched_assemblies:
                    unmatched_assembly_writer.writerow(
                        {"accession": accession, "organism_name": organism_name}
                    )
                if args.verbose:
                    print(
                        f"Warning: Organism name '{organism_name}' not found in taxonomy and no close matches found. Skipping assembly {accession}.",
                        file=sys.stderr,
                    )
                continue

        else:
            taxon_id = name_tax_id_lookup[organism_name]

        human_readable = tax_id_name_lookup[taxon_id]

        in_fasta = pyfastx.Fastx(assembly_path, uppercase=True, comment=True)

        for name, seq, comment in in_fasta:
            sequence_md5 = md5_checksum(seq)
            sequence_length = len(seq)

            metadata_record = {
                "unique_accession": name,
                "accession_description": comment,
                "taxon_id": taxon_id,
                "human_readable": human_readable,
                "source": "RefSeq",
                "source_accession": accession,
                "sequence_md5": sequence_md5,
                "sequence_length": sequence_length,
            }

            metadata_writer.writerow(metadata_record)

            fasta_writer.write(f">{name} {comment}\n{seq}\n")

    fasta_writer.close()


def main():
    parser = argparse.ArgumentParser(
        description="Create metadata TSV from genome assemblies and NCBI taxonomy."
    )
    parser.add_argument(
        "--taxonomy",
        type=str,
        required=True,
        help="Path to NCBI taxonomy directory containing names.dmp",
    )
    parser.add_argument(
        "--assembly-wildcard",
        type=str,
        required=True,
        help="Wildcard pattern to locate genome assembly directories.",
    )
    parser.add_argument(
        "--output_tsv", type=str, required=True, help="Output metadata TSV file."
    )
    parser.add_argument(
        "--output_fasta",
        type=str,
        required=True,
        help="Output combined FASTA file.",
    )
    parser.add_argument(
        "--unmatched_assemblies", type=str, help="File to log unmatched assemblies."
    )
    parser.add_argument(
        "-v", "--verbose", action="store_true", help="Enable verbose output."
    )
    args = parser.parse_args()
    run(args)


if __name__ == "__main__":
    main()
