import jsonschema
from pathlib import Path
import json
import csv
import pyfastx
import hashlib
import pandas as pd


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


def import_taxonomy_names(names_dmp_path: Path) -> dict:
    """Import taxonomy names from a names.dmp file.

    Args:
        names_dmp_path (Path): The path to the names.dmp file.

    Returns:
        dict: A dictionary mapping taxonomic IDs to their scientific names.
    """
    taxonomy = {}
    with open(names_dmp_path, "r") as f:
        for line in f:
            fields = [i.lstrip() for i in line.split("\t|")]
            taxon_id, name, name_type = fields[0], fields[1], fields[3]
            if "scientific name" in name_type:
                taxonomy[taxon_id] = name

    return taxonomy


def validate_metadata_tsv(metadata: Path, schema_path: Path) -> None:
    """Validate the metadata TSV file against the provided JSON schema.
    Args:
        metadata (Path): The path to the metadata TSV file.
        schema_path (Path): The path to the JSON schema file.

        Raises:
            jsonschema.ValidationError: If the data does not conform to the schema.
            jsonschema.SchemaError: If the schema itself is invalid.
    """
    with open(schema_path, "r") as f:
        schema = json.load(f)

    metadata_df = pd.read_csv(metadata, sep="\t")

    jsonschema.validate(instance=metadata_df.to_dict("records"), schema=schema)


def validate_index_json(index: Path, schema_path: Path) -> None:
    """Validate the database index JSON file against the provided JSON schema.
    Args:
        index (Path): The path to the database index JSON file.
        schema_path (Path): The path to the JSON schema file.

        Raises:
            jsonschema.ValidationError: If the data does not conform to the schema.
            jsonschema.SchemaError: If the schema itself is invalid.
    """
    with open(schema_path, "r") as f:
        schema = json.load(f)

    with open(index, "r") as f:
        index_data = json.load(f)

    jsonschema.validate(instance=index_data, schema=schema)


def run(args):

    name_lookup = import_taxonomy_names(args.names_dmp)

    with open(args.metadata, "r") as f:
        reader = csv.DictReader(f, delimiter="\t")

        metadata = [x for x in reader]

    # Basic structural validation of the metadata TSV
    validate_metadata_tsv(args.metadata, args.metadata_schema)
    print(f"Metadata TSV {args.metadata} passed schema validation.")

    # Basic structural validation of the database index JSON
    validate_index_json(args.index, args.index_schema)
    print(f"Database index JSON {args.index} passed schema validation.")

    # Check that all taxonomic IDs in the metadata are valid
    for row in metadata:
        taxon_id = row.get("taxon_id")
        if taxon_id not in name_lookup.keys():
            raise ValueError(
                f"taxon_id: {taxon_id} does not exist in the taxonomy database provided. Observed in record: {row}"
            )

        if row["human_readable"] != name_lookup[taxon_id]:
            raise ValueError(
                f"Mismatch between taxon_id: {taxon_id} and human_readable name: {row['human_readable']}."
            )

    # Validate FASTA file order, headers, and md5 checksums
    fasta = pyfastx.Fastx(str(args.fasta), uppercase=True, comment=True)

    fasta_entries = 0
    for metadata_row, (name, seq, comment) in zip(metadata, fasta):
        fasta_entries += 1
        if metadata_row["unique_accession"] != name:
            raise ValueError(
                f"FASTA header {name} does not match unique_accession {metadata_row['unique_accession']} in metadata. This may be due to ordering issues or incorrect metadata / FASTA header."
            )

        if metadata_row["accession_description"] != comment:
            raise ValueError(
                f"FASTA description for {name} does not match accession_description '{metadata_row['accession_description']}' in metadata."
            )

        if metadata_row["sequence_md5"] != md5_checksum(seq):
            raise ValueError(
                f"FASTA sequence for {name} does not match sequence_md5 '{metadata_row['sequence_md5']}' in metadata."
            )

        if int(metadata_row["sequence_length"]) != len(seq):
            raise ValueError(
                f"FASTA sequence length for {name}: {len(seq)} does not match sequence_length '{metadata_row['sequence_length']}' in metadata."
            )

    if fasta_entries != len(metadata):
        raise ValueError(
            f"Number of entries in FASTA ({fasta_entries}) does not match number of entries in metadata ({len(metadata)})."
        )

    print(f"Validation completed successfully. {fasta_entries} entries checked.")


def main():
    import argparse

    parser = argparse.ArgumentParser(
        description="Validate a metadata TSV file against a JSON schema and check taxonomic IDs."
    )
    parser.add_argument(
        "--names_dmp",
        type=Path,
        required=True,
        help="Path to the names.dmp file for taxonomy lookup.",
    )
    parser.add_argument(
        "--metadata_schema",
        type=Path,
        required=True,
        help="Path to the JSON schema file for metadata TSV validation.",
    )
    parser.add_argument(
        "--index_schema",
        type=Path,
        required=True,
        help="Path to the JSON schema file for database index JSON validation.",
    )
    parser.add_argument(
        "--metadata",
        type=Path,
        required=True,
        help="Path to the metadata TSV file to validate.",
    )
    parser.add_argument(
        "--index",
        type=Path,
        required=True,
        help="Path to the database index JSON file to validate.",
    )
    parser.add_argument(
        "--fasta",
        type=Path,
        required=True,
        help="Path to the FASTA file to validate.",
    )

    args = parser.parse_args()
    run(args)


if __name__ == "__main__":
    main()
