#!/usr/bin/env python3

"""Create safe Panaroo-generate-gffs inputs without altering source GFFs."""

import argparse
from pathlib import Path


def parse_arguments() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Copy reference GFFs into a new directory while changing only "
            "the legacy auxiliary attribute panaroo_ID= to panaroo_id=."
        )
    )
    parser.add_argument("input_dir", type=Path)
    parser.add_argument("output_dir", type=Path)
    return parser.parse_args()


def fasta_sequence(text: str, path: Path) -> tuple[str, str]:
    lines = text.splitlines()

    try:
        fasta_index = lines.index("##FASTA")
    except ValueError as error:
        raise ValueError(f"Missing ##FASTA directive: {path}") from error

    headers = [
        line[1:].split()[0]
        for line in lines[fasta_index + 1 :]
        if line.startswith(">")
    ]

    if len(headers) != 1:
        raise ValueError(
            f"Expected one FASTA sequence in {path}; found {len(headers)}"
        )

    sequence = "".join(
        line.strip()
        for line in lines[fasta_index + 1 :]
        if line and not line.startswith(">")
    ).upper()

    return headers[0], sequence


def main() -> None:
    args = parse_arguments()

    if not args.input_dir.is_dir():
        raise NotADirectoryError(args.input_dir)

    if args.output_dir.exists():
        raise FileExistsError(
            f"Refusing to overwrite existing output directory: "
            f"{args.output_dir}"
        )

    input_files = sorted(args.input_dir.glob("*.gff"))
    if not input_files:
        raise ValueError(f"No .gff files found in {args.input_dir}")

    prepared: list[tuple[Path, str]] = []
    fasta_ids: set[str] = set()
    replacements = 0

    for input_path in input_files:
        original = input_path.read_text(encoding="utf-8")
        fasta_id, original_sequence = fasta_sequence(original, input_path)

        if fasta_id in fasta_ids:
            raise ValueError(f"Duplicate FASTA identifier: {fasta_id}")
        fasta_ids.add(fasta_id)

        replacement_count = original.count("panaroo_ID=")
        converted = original.replace("panaroo_ID=", "panaroo_id=")
        replacements += replacement_count

        if "panaroo_ID=" in converted:
            raise RuntimeError(
                f"Uppercase panaroo_ID remains in {input_path}"
            )

        converted_id, converted_sequence = fasta_sequence(
            converted,
            input_path,
        )

        if fasta_id != converted_id:
            raise RuntimeError(f"FASTA identifier changed in {input_path}")

        if original_sequence != converted_sequence:
            raise RuntimeError(f"FASTA sequence changed in {input_path}")

        # Confirm that no field other than the exact auxiliary key changed.
        restored = converted.replace("panaroo_id=", "panaroo_ID=")
        if replacement_count and restored != original:
            raise RuntimeError(
                f"Unexpected content change while preparing {input_path}"
            )

        prepared.append((args.output_dir / input_path.name, converted))

    args.output_dir.mkdir(parents=True)

    for output_path, text in prepared:
        output_path.write_text(text, encoding="utf-8")

    output_files = sorted(args.output_dir.glob("*.gff"))
    if [path.name for path in output_files] != [
        path.name for path in input_files
    ]:
        raise RuntimeError("Output GFF filenames do not match input filenames")

    uppercase_remaining = sum(
        path.read_text(encoding="utf-8").count("panaroo_ID=")
        for path in output_files
    )
    if uppercase_remaining:
        raise RuntimeError(
            f"Uppercase panaroo_ID fields remain: {uppercase_remaining}"
        )

    print(f"Input GFFs: {len(input_files)}")
    print(f"Output GFFs: {len(output_files)}")
    print(f"panaroo_ID fields changed: {replacements}")
    print(f"FASTA identifiers: {len(fasta_ids)}")
    print(f"Wrote: {args.output_dir}")


if __name__ == "__main__":
    main()
