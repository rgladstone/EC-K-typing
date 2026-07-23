#!/usr/bin/env python3

"""Normalise one private Bakta-style GFF without modifying its sequence."""

import argparse
import re
from pathlib import Path


def parse_arguments() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Copy a single-sequence GFF into the EC-K-Typing naming scheme "
            "using a local-only numeric KL identifier."
        )
    )
    parser.add_argument("input_gff", type=Path)
    parser.add_argument("output_gff", type=Path)
    parser.add_argument("--group", required=True, choices=("2", "3"))
    parser.add_argument("--locus", required=True)
    parser.add_argument("--accession", default="PRIVATE")
    parser.add_argument(
        "--allow-atypical-boundaries",
        action="store_true",
        help=(
            "Accept boundaries that do not match the expected group 2 or "
            "group 3 orientation. Use only after manual inspection."
        ),
    )
    return parser.parse_args()


def validate_token(value: str, label: str) -> str:
    value = value.strip()
    if not re.fullmatch(r"[A-Za-z0-9.-]+", value):
        raise ValueError(
            f"{label} must contain only letters, numbers, dots or hyphens: "
            f"{value!r}"
        )
    return value


def parse_attributes(value: str) -> dict[str, str]:
    return {
        key: item
        for field in value.split(";")
        if "=" in field
        for key, item in [field.split("=", 1)]
    }


def boundary_text(attributes: dict[str, str]) -> str:
    keys = ("gene", "name", "Name", "product", "description")
    return " ".join(attributes.get(key, "") for key in keys).lower()


def matches_boundary(text: str, aliases: tuple[str, ...]) -> bool:
    return any(alias in text for alias in aliases)


def main() -> None:
    args = parse_arguments()

    if not args.input_gff.is_file():
        raise FileNotFoundError(args.input_gff)
    if args.output_gff.exists():
        raise FileExistsError(
            f"Refusing to overwrite existing output: {args.output_gff}"
        )

    locus = args.locus.strip().upper()
    if not re.fullmatch(r"KL\d+", locus):
        raise ValueError(
            "--locus must be a local numeric identifier such as KL9001"
        )

    accession = validate_token(args.accession, "accession")
    sequence_id = f"G{args.group}_{locus}_{accession}"
    lines = args.input_gff.read_text(encoding="utf-8").splitlines()

    try:
        fasta_index = lines.index("##FASTA")
    except ValueError as error:
        raise ValueError(
            f"Missing ##FASTA directive: {args.input_gff}"
        ) from error

    fasta_headers = [
        index
        for index, line in enumerate(lines[fasta_index + 1 :], fasta_index + 1)
        if line.startswith(">")
    ]
    if len(fasta_headers) != 1:
        raise ValueError(
            f"Expected one FASTA sequence; found {len(fasta_headers)}"
        )

    header_index = fasta_headers[0]
    sequence = "".join(
        line.strip()
        for line in lines[header_index + 1 :]
        if line and not line.startswith(">")
    ).upper()
    if not sequence:
        raise ValueError("The GFF FASTA sequence is empty")
    if "N" in sequence:
        raise ValueError("The private locus sequence contains N bases")
    if not re.fullmatch(r"[ACGT]+", sequence):
        raise ValueError("The private locus contains non-ACGT characters")

    output: list[str] = []
    cds_count = 0
    cds_boundaries: list[tuple[int, str]] = []

    for line_number, line in enumerate(lines[:fasta_index], 1):
        if line.startswith("##sequence-region "):
            output.append(f"##sequence-region {sequence_id} 1 {len(sequence)}")
            continue
        if not line or line.startswith("#"):
            output.append(line)
            continue

        fields = line.split("\t")
        if len(fields) != 9:
            raise ValueError(
                f"Expected 9 GFF fields at line {line_number}; "
                f"found {len(fields)}"
            )
        fields[0] = sequence_id

        try:
            start = int(fields[3])
            end = int(fields[4])
        except ValueError as error:
            raise ValueError(
                f"Invalid coordinates at line {line_number}"
            ) from error

        if start < 1 or end < start or end > len(sequence):
            raise ValueError(
                f"Feature outside the sequence at line {line_number}: "
                f"{start}-{end}"
            )
        if fields[2] == "CDS":
            cds_count += 1
            cds_boundaries.append(
                (
                    start,
                    boundary_text(parse_attributes(fields[8])),
                )
            )
        output.append("\t".join(fields))

    if cds_count == 0:
        raise ValueError("The private GFF contains no CDS features")

    cds_boundaries.sort()
    first_boundary = cds_boundaries[0][1]
    last_boundary = cds_boundaries[-1][1]
    if args.group == "2":
        first_aliases = (
            "kpsf",
            "arabinose 5-phosphate isomerase",
            "polysialic acid capsule expression protein",
        )
        last_aliases = (
            "kpsm",
            "polysialic acid transport protein",
            "abc transporter permease",
            "transport permease protein",
        )
        expected = "kpsF-to-kpsM"
    else:
        first_aliases = (
            "kpsd",
            "wza",
            "kpsm",
            "polysialic acid transport protein",
            "transport permease protein",
        )
        last_aliases = (
            "kpss",
            "capsule polysaccharide modification protein",
        )
        expected = "kpsD/kpsM-to-kpsS"

    boundaries_valid = (
        matches_boundary(first_boundary, first_aliases)
        and matches_boundary(last_boundary, last_aliases)
    )
    if not boundaries_valid and not args.allow_atypical_boundaries:
        raise ValueError(
            f"The CDS boundaries do not match the expected {expected} "
            f"group {args.group} orientation. First CDS: "
            f"{first_boundary!r}; last CDS: {last_boundary!r}. "
            "Inspect the locus and orientation, or use "
            "--allow-atypical-boundaries for a reviewed atypical locus."
        )

    output.extend(
        [
            "##FASTA",
            f">{sequence_id}",
            *(
                sequence[position : position + 60]
                for position in range(0, len(sequence), 60)
            ),
        ]
    )

    args.output_gff.parent.mkdir(parents=True, exist_ok=True)
    args.output_gff.write_text(
        "\n".join(output) + "\n",
        encoding="utf-8",
    )

    print(f"Input: {args.input_gff}")
    print(f"Output: {args.output_gff}")
    print(f"Sequence ID: {sequence_id}")
    print(f"Length: {len(sequence)}")
    print(f"CDS features: {cds_count}")
    print(
        "Boundary check: "
        + ("passed" if boundaries_valid else "explicitly overridden")
    )


if __name__ == "__main__":
    main()
