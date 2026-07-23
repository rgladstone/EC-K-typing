#!/usr/bin/env python3

"""Filter and stage Panaroo-generated GFFs using canonical reference names."""

import argparse
import re
from dataclasses import dataclass
from pathlib import Path
from typing import Optional


@dataclass
class ParsedGff:
    path: Path
    lines: list[str]
    feature_indexes: list[int]
    fasta_id: str
    sequence: str


def parse_arguments() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Remove reviewed short Panaroo refounds and newly assigned short "
            "hypothetical CDS annotations, validate the generated GFFs, and "
            "stage them using canonical reference filenames."
        )
    )
    parser.add_argument("generated_dir", type=Path)
    parser.add_argument("canonical_dir", type=Path)
    parser.add_argument("output_dir", type=Path)
    parser.add_argument(
        "--minimum-refound-length",
        type=int,
        default=200,
        help=(
            "Remove candidate_gene features shorter than this many bp "
            "(default: 200). Ordinary CDS features are never removed."
        ),
    )
    parser.add_argument(
        "--minimum-new-hypothetical-length",
        type=int,
        default=200,
        help=(
            "Remove newly assigned hypothetical CDS features shorter than "
            "this many bp (default: 200)."
        ),
    )
    parser.add_argument(
        "--first-new-db-group",
        type=int,
        default=257,
        help=(
            "First DB_group number assigned during this database expansion "
            "(default: 257). Lower-numbered historical groups are retained."
        ),
    )
    parser.add_argument(
        "--remove-private-short-hypotheticals",
        action="store_true",
        help=(
            "Also remove hypothetical CDS features shorter than "
            "--minimum-new-hypothetical-length when their stable name "
            "starts with private_group_. Intended only for private builds."
        ),
    )
    return parser.parse_args()


def parse_attributes(value: str, path: Path, line_number: int) -> dict[str, str]:
    attributes: dict[str, str] = {}

    for field in value.split(";"):
        if not field:
            continue
        if "=" not in field:
            raise ValueError(
                f"Malformed attribute in {path}:{line_number}: {field!r}"
            )
        key, item = field.split("=", 1)
        if not key:
            raise ValueError(
                f"Empty attribute key in {path}:{line_number}"
            )
        attributes[key] = item

    return attributes


def parse_gff(path: Path) -> ParsedGff:
    lines = path.read_text(encoding="utf-8").splitlines()

    try:
        fasta_index = lines.index("##FASTA")
    except ValueError as error:
        raise ValueError(f"Missing ##FASTA directive: {path}") from error

    feature_indexes: list[int] = []

    for index, line in enumerate(lines[:fasta_index]):
        if not line or line.startswith("#"):
            continue

        fields = line.split("\t")
        if len(fields) != 9:
            raise ValueError(
                f"Expected 9 GFF columns in {path}:{index + 1}; "
                f"found {len(fields)}"
            )

        try:
            start = int(fields[3])
            end = int(fields[4])
        except ValueError as error:
            raise ValueError(
                f"Invalid feature coordinates in {path}:{index + 1}"
            ) from error

        if start < 1 or end < start:
            raise ValueError(
                f"Invalid feature interval in {path}:{index + 1}: "
                f"{start}-{end}"
            )

        if fields[6] not in {"+", "-"}:
            raise ValueError(
                f"Invalid feature strand in {path}:{index + 1}: {fields[6]!r}"
            )

        parse_attributes(fields[8], path, index + 1)
        feature_indexes.append(index)

    fasta_lines = lines[fasta_index + 1 :]
    headers = [
        (index, line)
        for index, line in enumerate(fasta_lines)
        if line.startswith(">")
    ]

    if len(headers) != 1:
        raise ValueError(
            f"Expected one FASTA sequence in {path}; found {len(headers)}"
        )

    header_index, header = headers[0]
    fasta_id = header[1:].split()[0]

    if not fasta_id:
        raise ValueError(f"Empty FASTA identifier in {path}")

    sequence = "".join(
        line.strip()
        for line in fasta_lines[header_index + 1 :]
        if line and not line.startswith(">")
    ).upper()

    if not sequence:
        raise ValueError(f"Empty FASTA sequence in {path}")

    if "N" in sequence:
        raise ValueError(f"FASTA sequence contains N bases: {path}")

    for index in feature_indexes:
        fields = lines[index].split("\t")
        if int(fields[4]) > len(sequence):
            raise ValueError(
                f"Feature extends beyond the sequence in {path}:{index + 1}"
            )

    return ParsedGff(
        path=path,
        lines=lines,
        feature_indexes=feature_indexes,
        fasta_id=fasta_id,
        sequence=sequence,
    )


def load_directory(path: Path) -> dict[str, ParsedGff]:
    if not path.is_dir():
        raise NotADirectoryError(path)

    files = sorted(path.glob("*.gff"))
    if not files:
        raise ValueError(f"No .gff files found in {path}")

    records: dict[str, ParsedGff] = {}

    for file_path in files:
        record = parse_gff(file_path)
        if record.fasta_id in records:
            raise ValueError(
                f"Duplicate FASTA identifier {record.fasta_id!r}: "
                f"{records[record.fasta_id].path} and {file_path}"
            )
        records[record.fasta_id] = record

    return records


def is_is_related(attributes: dict[str, str]) -> bool:
    name = attributes.get("name", attributes.get("Name", "")).lower()
    description = attributes.get(
        "description", attributes.get("product", "")
    ).lower()
    dbxref = attributes.get("Dbxref", attributes.get("db_xref", ""))

    return (
        name == "tnp"
        or "transposase" in description
        or "insertion sequence" in description
        or bool(re.search(r"(^|,)IS:", dbxref))
    )


def db_group_number(name: str) -> Optional[int]:
    match = re.fullmatch(r"DB_group_(\d+)", name)
    return int(match.group(1)) if match else None


def filtered_lines(
    record: ParsedGff,
    minimum_refound_length: int,
    minimum_new_hypothetical_length: int,
    first_new_db_group: int,
    remove_private_short_hypotheticals: bool,
) -> tuple[
    list[str],
    list[tuple[str, int, str]],
    list[tuple[str, int, str]],
    list[tuple[str, int, str]],
]:
    remove_indexes: set[int] = set()
    removed_refounds: list[tuple[str, int, str]] = []
    removed_new_hypotheticals: list[tuple[str, int, str]] = []
    retained_candidates: list[tuple[str, int, str]] = []

    for index in record.feature_indexes:
        fields = record.lines[index].split("\t")
        feature_type = fields[2]
        length = int(fields[4]) - int(fields[3]) + 1
        attributes = parse_attributes(fields[8], record.path, index + 1)
        name = attributes.get("name", attributes.get("Name", ""))
        description = attributes.get(
            "description", attributes.get("product", "")
        )

        if is_is_related(attributes):
            raise ValueError(
                f"IS-related feature remains in {record.path}:{index + 1}"
            )

        if "~~~" in name or "~~~" in attributes.get(
            "panaroo_gene_cluster", ""
        ):
            raise ValueError(
                f"Merged cluster name remains in {record.path}:{index + 1}"
            )

        if name in {"", "No_name"}:
            raise ValueError(
                f"Missing curated gene name in {record.path}:{index + 1}"
            )

        if feature_type == "candidate_gene":
            if length < minimum_refound_length:
                remove_indexes.add(index)
                removed_refounds.append((record.fasta_id, length, name))
            else:
                retained_candidates.append((record.fasta_id, length, name))
            continue

        group_number = db_group_number(name)
        is_private_group = bool(
            re.fullmatch(
                r"private_(?:group|composite|split)_[0-9a-f]{12}",
                name,
            )
        )
        if (
            feature_type == "CDS"
            and length < minimum_new_hypothetical_length
            and description.strip().lower() == "hypothetical protein"
            and (
                (
                    group_number is not None
                    and group_number >= first_new_db_group
                )
                or (
                    remove_private_short_hypotheticals
                    and is_private_group
                )
            )
        ):
            remove_indexes.add(index)
            removed_new_hypotheticals.append(
                (record.fasta_id, length, name)
            )

    output_lines = [
        line
        for index, line in enumerate(record.lines)
        if index not in remove_indexes
    ]

    return (
        output_lines,
        removed_refounds,
        removed_new_hypotheticals,
        retained_candidates,
    )


def main() -> None:
    args = parse_arguments()

    if args.minimum_refound_length < 1:
        raise ValueError("--minimum-refound-length must be positive")
    if args.minimum_new_hypothetical_length < 1:
        raise ValueError(
            "--minimum-new-hypothetical-length must be positive"
        )
    if args.first_new_db_group < 1:
        raise ValueError("--first-new-db-group must be positive")

    if args.output_dir.exists():
        raise FileExistsError(
            f"Refusing to overwrite existing output directory: "
            f"{args.output_dir}"
        )

    generated = load_directory(args.generated_dir)
    canonical = load_directory(args.canonical_dir)

    missing = set(canonical) - set(generated)
    extra = set(generated) - set(canonical)

    if missing:
        raise ValueError(
            "Generated GFFs are missing FASTA identifiers: "
            + ", ".join(sorted(missing))
        )

    if extra:
        raise ValueError(
            "Generated GFFs contain unexpected FASTA identifiers: "
            + ", ".join(sorted(extra))
        )

    prepared: list[tuple[Path, list[str]]] = []
    removed_refounds: list[tuple[str, int, str]] = []
    removed_new_hypotheticals: list[tuple[str, int, str]] = []
    retained_candidates: list[tuple[str, int, str]] = []
    output_names: set[str] = set()

    for fasta_id in sorted(generated):
        generated_record = generated[fasta_id]
        canonical_record = canonical[fasta_id]

        if generated_record.sequence != canonical_record.sequence:
            raise ValueError(
                f"FASTA sequence changed for {fasta_id}: "
                f"{canonical_record.path} versus {generated_record.path}"
            )

        output_name = canonical_record.path.name
        if output_name in output_names:
            raise ValueError(f"Duplicate output filename: {output_name}")
        output_names.add(output_name)

        (
            lines,
            removed_refounds_here,
            removed_hypotheticals_here,
            retained_here,
        ) = filtered_lines(
            generated_record,
            args.minimum_refound_length,
            args.minimum_new_hypothetical_length,
            args.first_new_db_group,
            args.remove_private_short_hypotheticals,
        )
        prepared.append((args.output_dir / output_name, lines))
        removed_refounds.extend(removed_refounds_here)
        removed_new_hypotheticals.extend(
            removed_hypotheticals_here
        )
        retained_candidates.extend(retained_here)

    args.output_dir.mkdir(parents=True)

    for output_path, lines in prepared:
        output_path.write_text(
            "\n".join(lines) + "\n",
            encoding="utf-8",
        )

    staged = load_directory(args.output_dir)

    if set(staged) != set(canonical):
        raise RuntimeError("Staged FASTA identifiers changed after writing")

    for fasta_id in staged:
        if staged[fasta_id].sequence != canonical[fasta_id].sequence:
            raise RuntimeError(
                f"Staged FASTA sequence changed for {fasta_id}"
            )

    remaining_short_candidates = []
    for record in staged.values():
        for index in record.feature_indexes:
            fields = record.lines[index].split("\t")
            if fields[2] != "candidate_gene":
                continue
            length = int(fields[4]) - int(fields[3]) + 1
            if length < args.minimum_refound_length:
                remaining_short_candidates.append(
                    (record.fasta_id, length)
                )

    if remaining_short_candidates:
        raise RuntimeError(
            f"Short candidate genes remain: {remaining_short_candidates}"
        )

    remaining_new_short_hypotheticals = []
    for record in staged.values():
        for index in record.feature_indexes:
            fields = record.lines[index].split("\t")
            if fields[2] != "CDS":
                continue
            length = int(fields[4]) - int(fields[3]) + 1
            attributes = parse_attributes(
                fields[8], record.path, index + 1
            )
            name = attributes.get("name", attributes.get("Name", ""))
            description = attributes.get(
                "description", attributes.get("product", "")
            )
            group_number = db_group_number(name)
            is_private_group = bool(
                re.fullmatch(
                    r"private_(?:group|composite|split)_[0-9a-f]{12}",
                    name,
                )
            )
            if (
                length < args.minimum_new_hypothetical_length
                and description.strip().lower() == "hypothetical protein"
                and (
                    (
                        group_number is not None
                        and group_number >= args.first_new_db_group
                    )
                    or (
                        args.remove_private_short_hypotheticals
                        and is_private_group
                    )
                )
            ):
                remaining_new_short_hypotheticals.append(
                    (record.fasta_id, length, name)
                )

    if remaining_new_short_hypotheticals:
        raise RuntimeError(
            "New short hypothetical CDS features remain: "
            f"{remaining_new_short_hypotheticals}"
        )

    print(f"Generated GFFs: {len(generated)}")
    print(f"Canonical GFFs: {len(canonical)}")
    print(f"Staged GFFs: {len(staged)}")
    print(
        f"Short candidate genes removed: {len(removed_refounds)}"
    )
    print(
        "New short hypothetical CDS features removed: "
        f"{len(removed_new_hypotheticals)}"
    )
    print(f"Candidate genes retained: {len(retained_candidates)}")

    for fasta_id, length, name in retained_candidates:
        print(
            f"RETAINED\t{fasta_id}\t{length}\t{name}"
        )

    print(f"Wrote: {args.output_dir}")


if __name__ == "__main__":
    main()
