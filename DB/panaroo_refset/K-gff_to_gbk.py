#!/usr/bin/env python3

"""Convert curated EC-K-Typing GFFs into a Kaptive GenBank database."""

import argparse
import csv
import datetime
import re
from dataclasses import dataclass
from io import StringIO
from pathlib import Path

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqFeature import FeatureLocation, SeqFeature
from Bio.SeqRecord import SeqRecord


METADATA_COLUMNS = {"locus", "phenotype", "note"}


@dataclass
class LocusMetadata:
    phenotype: str = ""
    note: str = ""


@dataclass
class ParsedGff:
    path: Path
    locus: str
    phenotype: str
    accession: str
    sequence: str
    features: list[SeqFeature]
    skipped_candidates: list[tuple[str, int]]


def parse_arguments() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Convert curated Panaroo GFF records to a multi-record GenBank "
            "database compatible with Kaptive."
        )
    )
    parser.add_argument(
        "gffs",
        nargs="+",
        type=Path,
        help="Curated GFF files to include in the database.",
    )
    parser.add_argument(
        "-o",
        "--output",
        required=True,
        type=Path,
        help="Output GenBank database path.",
    )
    parser.add_argument(
        "--metadata",
        type=Path,
        help=(
            "Optional CSV containing exactly: locus, phenotype, note. "
            "Phenotype values override filename-derived K types."
        ),
    )
    parser.add_argument(
        "--locus-version",
        default="v2",
        help="Per-locus VERSION suffix (default: v2).",
    )
    parser.add_argument(
        "--line-ending",
        choices=("crlf", "lf"),
        default="crlf",
        help=(
            "Output line ending. CRLF preserves the existing repository "
            "database format (default: crlf)."
        ),
    )
    return parser.parse_args()


def sanitize(value: str) -> str:
    value = re.sub(r"\s+", " ", str(value)).strip()
    return "".join(
        character
        for character in value
        if 32 <= ord(character) <= 126
    )


def parse_attributes(
    value: str,
    path: Path,
    line_number: int,
) -> dict[str, str]:
    attributes: dict[str, str] = {}

    for field in value.split(";"):
        if not field:
            continue
        if "=" not in field:
            raise ValueError(
                f"Malformed attribute in {path}:{line_number}: {field!r}"
            )
        key, item = field.split("=", 1)
        attributes[key] = item

    return attributes


def filename_metadata(path: Path) -> tuple[str, str]:
    locus_match = re.search(r"(KL\d+)", path.name, re.IGNORECASE)
    if not locus_match:
        raise ValueError(f"Cannot derive KL identifier from {path.name}")

    locus = locus_match.group(1).upper()
    phenotype_match = re.match(
        r"(K[A-Z0-9-]+)_G[23]_",
        path.name,
        re.IGNORECASE,
    )
    phenotype = phenotype_match.group(1).upper() if phenotype_match else ""

    return locus, phenotype


def accession_from_fasta_id(fasta_id: str, locus: str, path: Path) -> str:
    match = re.search(
        rf"(?:^|_){re.escape(locus)}_(.+)$",
        fasta_id,
        re.IGNORECASE,
    )

    if not match:
        raise ValueError(
            f"Cannot derive accession after {locus} from FASTA identifier "
            f"{fasta_id!r} in {path}"
        )

    accession = sanitize(match.group(1))
    if not accession:
        raise ValueError(f"Empty accession in FASTA identifier: {fasta_id}")

    return accession


def parse_gff(path: Path) -> ParsedGff:
    if not path.is_file():
        raise FileNotFoundError(path)

    locus, phenotype = filename_metadata(path)
    features: list[SeqFeature] = []
    fasta_headers: list[str] = []
    sequence_parts: list[str] = []
    skipped_candidates: list[tuple[str, int]] = []
    in_fasta = False

    with path.open(encoding="utf-8") as handle:
        for line_number, raw_line in enumerate(handle, 1):
            line = raw_line.rstrip("\r\n")

            if line == "##FASTA":
                in_fasta = True
                continue

            if in_fasta:
                if line.startswith(">"):
                    fasta_headers.append(line[1:].split()[0])
                elif line:
                    sequence_parts.append(line.strip())
                continue

            if not line or line.startswith("#"):
                continue

            fields = line.split("\t")
            if len(fields) != 9:
                raise ValueError(
                    f"Expected 9 GFF columns in {path}:{line_number}; "
                    f"found {len(fields)}"
                )

            feature_type = fields[2]
            if feature_type not in {"CDS", "candidate_gene"}:
                continue

            try:
                start = int(fields[3])
                end = int(fields[4])
            except ValueError as error:
                raise ValueError(
                    f"Invalid coordinates in {path}:{line_number}"
                ) from error

            if start < 1 or end < start:
                raise ValueError(
                    f"Invalid interval in {path}:{line_number}: "
                    f"{start}-{end}"
                )

            strand_text = fields[6]
            if strand_text not in {"+", "-"}:
                raise ValueError(
                    f"Invalid strand in {path}:{line_number}: "
                    f"{strand_text!r}"
                )

            attributes = parse_attributes(
                fields[8],
                path,
                line_number,
            )
            gene = sanitize(
                attributes.get("name", attributes.get("Name", ""))
            )
            product = sanitize(
                attributes.get(
                    "description",
                    attributes.get("product", ""),
                )
            )

            if gene in {"", "No_name"}:
                raise ValueError(
                    f"Missing curated gene name in {path}:{line_number}"
                )

            if "~~~" in gene:
                raise ValueError(
                    f"Merged gene name remains in {path}:{line_number}"
                )

            feature_length = end - start + 1
            if feature_length % 3:
                if feature_type == "candidate_gene":
                    skipped_candidates.append((gene, feature_length))
                    continue
                raise ValueError(
                    f"CDS length is not a multiple of three in "
                    f"{path}:{line_number}: {gene} ({feature_length} bp)"
                )

            qualifiers = {
                "gene": gene,
                "product": product or "hypothetical protein",
            }

            features.append(
                SeqFeature(
                    FeatureLocation(
                        start - 1,
                        end,
                        strand=1 if strand_text == "+" else -1,
                    ),
                    type="CDS",
                    qualifiers=qualifiers,
                )
            )

    if len(fasta_headers) != 1:
        raise ValueError(
            f"Expected one FASTA sequence in {path}; "
            f"found {len(fasta_headers)}"
        )

    sequence = "".join(sequence_parts).upper()
    if not sequence:
        raise ValueError(f"Empty FASTA sequence in {path}")

    if "N" in sequence:
        raise ValueError(f"FASTA sequence contains N bases: {path}")

    if any(int(feature.location.end) > len(sequence) for feature in features):
        raise ValueError(f"A feature extends beyond the FASTA sequence: {path}")

    fasta_id = fasta_headers[0]
    fasta_locus_match = re.search(r"(KL\d+)", fasta_id, re.IGNORECASE)

    if not fasta_locus_match:
        raise ValueError(
            f"Cannot derive KL identifier from FASTA ID {fasta_id!r}"
        )

    fasta_locus = fasta_locus_match.group(1).upper()
    if fasta_locus != locus:
        raise ValueError(
            f"Filename locus {locus} does not match FASTA locus "
            f"{fasta_locus} in {path}"
        )

    accession = accession_from_fasta_id(fasta_id, locus, path)

    return ParsedGff(
        path=path,
        locus=locus,
        phenotype=phenotype,
        accession=accession,
        sequence=sequence,
        features=features,
        skipped_candidates=skipped_candidates,
    )


def load_metadata(path: Path | None) -> dict[str, LocusMetadata]:
    if path is None:
        return {}

    if not path.is_file():
        raise FileNotFoundError(path)

    with path.open(newline="", encoding="utf-8-sig") as handle:
        reader = csv.DictReader(handle)
        if set(reader.fieldnames or []) != METADATA_COLUMNS:
            raise ValueError(
                f"{path} must contain exactly these columns: "
                + ", ".join(sorted(METADATA_COLUMNS))
            )
        rows = list(reader)

    metadata: dict[str, LocusMetadata] = {}

    for row in rows:
        locus = row["locus"].strip().upper()
        if not re.fullmatch(r"KL\d+", locus):
            raise ValueError(f"Invalid locus in {path}: {locus!r}")
        if locus in metadata:
            raise ValueError(f"Duplicate metadata locus: {locus}")

        metadata[locus] = LocusMetadata(
            phenotype=sanitize(row["phenotype"]),
            note=sanitize(row["note"]),
        )

    return metadata


def add_locus_tags(locus: ParsedGff) -> None:
    for position, feature in enumerate(locus.features, 1):
        feature.qualifiers["locus_tag"] = (
            f"{locus.locus}_{position:05d}"
        )


def make_record(
    locus: ParsedGff,
    metadata: LocusMetadata,
) -> SeqRecord:
    phenotype = metadata.phenotype or locus.phenotype
    description = (
        "Escherichia coli DNA, capsular polysaccharide synthesis "
        "gene cluster"
    )
    if phenotype:
        description += f", serotype: {phenotype}"

    record = SeqRecord(
        Seq(locus.sequence),
        id=locus.accession,
        name=f"{locus.locus}_{locus.accession}",
        description=description,
    )
    record.annotations["molecule_type"] = "DNA"
    record.annotations["source"] = "Escherichia coli"
    record.annotations["organism"] = "Escherichia coli"
    record.annotations["taxonomy"] = [
        "Bacteria",
        "Proteobacteria",
        "Gammaproteobacteria",
        "Enterobacteriales",
        "Enterobacteriaceae",
        "Escherichia",
    ]

    notes = [f"K locus:{locus.locus}"]
    if phenotype:
        # Kaptive's default type regex captures punctuation and spaces only
        # when the value is separated from "type:" by a space.
        notes.append(f"K type: {phenotype}")
    if metadata.note:
        notes.append(metadata.note)
    notes.append("IS-element related annotations removed")

    source_qualifiers: dict[str, object] = {
        "organism": "Escherichia coli",
        "mol_type": "genomic DNA",
        "db_xref": "taxon:562",
        "note": notes,
    }

    if phenotype.upper().startswith("K"):
        source_qualifiers["serotype"] = phenotype

    source = SeqFeature(
        FeatureLocation(0, len(locus.sequence), strand=1),
        type="source",
        qualifiers=source_qualifiers,
    )

    add_locus_tags(locus)
    record.features = [source, *locus.features]
    return record


def format_record(
    record: SeqRecord,
    locus: str,
    accession: str,
    locus_version: str,
    date: str,
) -> str:
    buffer = StringIO()
    SeqIO.write(record, buffer, "genbank")
    lines = buffer.getvalue().splitlines()
    output: list[str] = []
    version_written = False
    comment_written = False
    comment = [
        "COMMENT     Annotation reconciled using Bakta v1.10.x, Panaroo",
        "            v1.5.2 for clustering and Panaroo v1.6.0 for GFF",
        "            generation.",
        "            https://github.com/oschwengers/bakta",
        "            https://github.com/gtonkinhill/panaroo",
    ]

    for line in lines:
        if line.startswith("LOCUS"):
            locus_name = f"{locus}_{accession}"
            padding = max(1, 45 - (9 + len(locus_name)))
            length = len(record.seq)
            output.append(
                f"LOCUS       {locus_name}"
                f"{' ' * padding}{length:>6} bp    DNA     linear   "
                f"UNK {date}"
            )
        elif line.startswith("ACCESSION"):
            output.append(line)
            output.append(f"VERSION     {locus}_{locus_version}")
            version_written = True
        elif line.startswith("VERSION") and version_written:
            continue
        elif (
            line.startswith(("KEYWORDS", "SOURCE", "FEATURES"))
            and not comment_written
        ):
            output.extend(comment)
            output.append(line)
            comment_written = True
        elif line.startswith("//"):
            output.append("//")
        else:
            output.append(line)

    return "\n".join(output)


def main() -> None:
    args = parse_arguments()

    if args.output.exists():
        raise FileExistsError(
            f"Refusing to overwrite existing output: {args.output}"
        )

    metadata = load_metadata(args.metadata)
    loci = [parse_gff(path) for path in args.gffs]

    locus_names = [locus.locus for locus in loci]
    accessions = [locus.accession for locus in loci]

    if len(locus_names) != len(set(locus_names)):
        duplicates = sorted(
            name
            for name in set(locus_names)
            if locus_names.count(name) > 1
        )
        raise ValueError(
            "Duplicate KL identifiers: " + ", ".join(duplicates)
        )

    if len(accessions) != len(set(accessions)):
        duplicates = sorted(
            accession
            for accession in set(accessions)
            if accessions.count(accession) > 1
        )
        raise ValueError(
            "Duplicate accessions: " + ", ".join(duplicates)
        )

    unknown_metadata = set(metadata) - set(locus_names)
    if unknown_metadata:
        raise ValueError(
            "Metadata loci absent from GFF inputs: "
            + ", ".join(sorted(unknown_metadata))
        )

    date = datetime.date.today().strftime("%d-%b-%Y").upper()
    records = [
        format_record(
            make_record(locus, metadata.get(locus.locus, LocusMetadata())),
            locus.locus,
            locus.accession,
            args.locus_version,
            date,
        )
        for locus in loci
    ]

    text = "\n".join(records)
    if args.line_ending == "crlf":
        text = text.replace("\n", "\r\n")

    args.output.parent.mkdir(parents=True, exist_ok=True)
    with args.output.open("w", encoding="ascii", newline="") as handle:
        handle.write(text)

    print(f"Records: {len(loci)}")
    print(f"CDS features: {sum(len(locus.features) for locus in loci)}")
    print(
        "Invalid candidate features skipped: "
        + str(sum(len(locus.skipped_candidates) for locus in loci))
    )
    for locus in loci:
        for gene, length in locus.skipped_candidates:
            print(
                f"SKIPPED\t{locus.locus}\t{gene}\t{length}\t"
                "not a multiple of three"
            )
    print(
        "Phenotype records: "
        + str(
            sum(
                bool(
                    metadata.get(locus.locus, LocusMetadata()).phenotype
                    or locus.phenotype
                )
                for locus in loci
            )
        )
    )
    print(f"Wrote: {args.output}")


if __name__ == "__main__":
    main()
