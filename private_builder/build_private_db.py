#!/usr/bin/env python3

"""Build and self-test a local EC-K-Typing database with one private locus."""

import argparse
import csv
import hashlib
import json
import os
import re
import shlex
import shutil
import subprocess
import sys
from datetime import datetime, timezone
from pathlib import Path


BUILDER_VERSION = "1.0.7"


def parse_arguments() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Add one private Bakta GFF to the frozen EC-K-Typing reference "
            "set and produce a validated local Kaptive database. The source "
            "repository is never modified."
        )
    )
    parser.add_argument("--input-gff", required=True, type=Path)
    parser.add_argument("--output-dir", required=True, type=Path)
    parser.add_argument("--group", required=True, choices=("2", "3"))
    parser.add_argument("--locus", required=True)
    parser.add_argument("--accession", default="PRIVATE")
    parser.add_argument("--phenotype", default="")
    parser.add_argument("--note", default="")
    parser.add_argument(
        "--allow-atypical-boundaries",
        action="store_true",
        help=(
            "Accept non-standard group boundaries after manual inspection."
        ),
    )
    parser.add_argument(
        "--threads",
        type=int,
        default=min(4, os.cpu_count() or 1),
    )
    parser.add_argument(
        "--reference-root",
        type=Path,
        default=Path(
            os.environ.get("EC_K_REFERENCE_ROOT", "/opt/ec-k-typing")
        ),
    )
    parser.add_argument(
        "--keep-work",
        action="store_true",
        help="Keep intermediate Panaroo outputs after a successful build.",
    )
    return parser.parse_args()


def environment_prefix(variable: str, environment: str) -> list[str]:
    override = os.environ.get(variable, "").strip()
    if override:
        return shlex.split(override)

    root = Path(os.environ.get("MAMBA_ROOT_PREFIX", "/opt/conda"))
    environment_bin = root / "envs" / environment / "bin"
    if not environment_bin.is_dir():
        raise NotADirectoryError(environment_bin)

    path = f"{environment_bin}:{os.environ.get('PATH', '')}"
    return ["/usr/bin/env", f"PATH={path}"]


def run(command: list[object], cwd: Path) -> None:
    printable = " ".join(shlex.quote(str(item)) for item in command)
    print(f"+ {printable}", flush=True)
    subprocess.run(
        [str(item) for item in command],
        cwd=cwd,
        check=True,
    )


def capture(command: list[object], cwd: Path) -> str:
    completed = subprocess.run(
        [str(item) for item in command],
        cwd=cwd,
        check=True,
        text=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
    )
    return " ".join(completed.stdout.split())


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def reference_sha256(paths: list[Path], root: Path) -> str:
    digest = hashlib.sha256()
    for path in sorted(paths, key=lambda item: str(item.relative_to(root))):
        relative = str(path.relative_to(root)).encode("utf-8")
        digest.update(relative)
        digest.update(b"\0")
        with path.open("rb") as handle:
            for chunk in iter(lambda: handle.read(1024 * 1024), b""):
                digest.update(chunk)
        digest.update(b"\0")
    return digest.hexdigest()


def is_within(path: Path, directory: Path) -> bool:
    try:
        path.relative_to(directory)
        return True
    except ValueError:
        return False


def clean_text(value: str, label: str, maximum: int = 200) -> str:
    value = value.strip()
    if any(ord(character) < 32 for character in value):
        raise ValueError(f"{label} must not contain control characters")
    if len(value) > maximum:
        raise ValueError(f"{label} must be at most {maximum} characters")
    return value


def reference_version(path: Path) -> str:
    text = path.read_text(encoding="utf-8")
    match = re.search(r'^version\s*=\s*"([^"]+)"', text, re.MULTILINE)
    if not match:
        raise ValueError(f"Cannot read version from {path}")
    return match.group(1)


def private_metadata(
    base_path: Path,
    output_path: Path,
    locus: str,
    phenotype: str,
    note: str,
) -> None:
    with base_path.open(newline="", encoding="utf-8-sig") as handle:
        reader = csv.DictReader(handle)
        rows = list(reader)
        fieldnames = reader.fieldnames

    if fieldnames != ["locus", "phenotype", "note"]:
        raise ValueError(
            f"Unexpected metadata columns in {base_path}: {fieldnames}"
        )
    if any(row["locus"].upper() == locus for row in rows):
        raise ValueError(f"Private locus already exists in metadata: {locus}")

    rows.append(
        {
            "locus": locus,
            "phenotype": phenotype,
            "note": note,
        }
    )
    with output_path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)


def validate_self_test(
    path: Path,
    expected_rows: int,
    private_locus: str,
    expected_phenotype: str,
) -> dict[str, int]:
    with path.open(newline="", encoding="utf-8") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))

    problems: list[str] = []
    if len(rows) != expected_rows:
        problems.append(
            f"Expected {expected_rows} result rows; found {len(rows)}"
        )

    assemblies = [row["Assembly"] for row in rows]
    if len(assemblies) != len(set(assemblies)):
        problems.append("Self-test contains duplicate assembly names")

    for row in rows:
        assembly = row["Assembly"]
        if row["Best match locus"] != assembly:
            problems.append(
                f"{assembly}: best locus is {row['Best match locus']}"
            )
        if row["Identity"] != "100.00%":
            problems.append(
                f"{assembly}: identity is {row['Identity']}"
            )
        if row["Coverage"] != "100.00%":
            problems.append(
                f"{assembly}: coverage is {row['Coverage']}"
            )
        if row["Match confidence"] == "Untypeable":
            problems.append(f"{assembly}: self-reference is untypeable")

    private_rows = [
        row for row in rows if row["Assembly"] == private_locus
    ]
    if len(private_rows) != 1:
        problems.append(
            f"Expected one self-test row for {private_locus}; "
            f"found {len(private_rows)}"
        )
    elif (
        expected_phenotype
        and private_rows[0]["Best match type"] != expected_phenotype
    ):
        problems.append(
            f"{private_locus}: phenotype is "
            f"{private_rows[0]['Best match type']!r}, expected "
            f"{expected_phenotype!r}"
        )

    if problems:
        raise RuntimeError(
            "Private database self-test failed:\n"
            + "\n".join(f"- {problem}" for problem in problems)
        )

    return {
        "rows": len(rows),
        "problem_flags": sum(bool(row["Problems"]) for row in rows),
        "untypeable": sum(
            row["Match confidence"] == "Untypeable" for row in rows
        ),
    }


def main() -> None:
    args = parse_arguments()
    root = args.reference_root.resolve()
    input_gff = args.input_gff.resolve()
    output_dir = args.output_dir.resolve()

    if not input_gff.is_file():
        raise FileNotFoundError(input_gff)
    if output_dir.exists():
        raise FileExistsError(
            f"Refusing to overwrite existing output directory: {output_dir}"
        )
    if is_within(output_dir, root):
        raise ValueError(
            "The output directory must be outside the frozen reference "
            "repository"
        )
    if args.threads < 1:
        raise ValueError("--threads must be positive")

    locus = args.locus.strip().upper()
    if not re.fullmatch(r"KL\d+", locus):
        raise ValueError(
            "--locus must be a local numeric identifier such as KL9001"
        )
    if int(locus[2:]) < 9000:
        raise ValueError(
            "Private locus numbers must be KL9000 or greater to avoid "
            "confusion with official EC-K-Typing identifiers"
        )
    accession = clean_text(args.accession, "accession", maximum=80)
    if not re.fullmatch(r"[A-Za-z0-9.-]+", accession):
        raise ValueError(
            "accession must contain only letters, numbers, dots or hyphens"
        )
    phenotype = clean_text(args.phenotype, "phenotype")
    note = clean_text(args.note, "note", maximum=500)

    reference_gffs = root / "DB" / "panaroo_refset" / "gffs"
    scripts = root / "DB" / "panaroo_refset"
    base_metadata = scripts / "locus_metadata.csv"
    config = root / "EC-K-typing_group2and3.toml"

    required = [
        reference_gffs,
        scripts / "edit_novel_gff.sh",
        scripts / "prepare_gffs_for_generate.py",
        scripts / "finalise_generated_gffs.py",
        scripts / "K-gff_to_gbk.py",
        base_metadata,
        config,
        root / "private_builder" / "normalise_private_gff.py",
        root / "private_builder" / "reconcile_graph.py",
    ]
    for path in required:
        if not path.exists():
            raise FileNotFoundError(path)

    official_gffs = sorted(reference_gffs.glob("*.gff"))
    if not official_gffs:
        raise ValueError(f"No official reference GFFs in {reference_gffs}")
    if any(
        re.search(rf"(?:^|_){locus}(?:_|$)", path.name)
        for path in official_gffs
    ):
        raise ValueError(f"Private locus collides with an official locus: {locus}")

    output_dir.mkdir(parents=True)
    work = output_dir / ".work"
    work.mkdir()
    inputs = work / "inputs"
    inputs.mkdir()

    builder_run = environment_prefix(
        "EC_K_BUILDER_RUN",
        "builder",
    )
    panaroo_152 = environment_prefix(
        "EC_K_PANAROO_152_RUN",
        "panaroo_1.5.2",
    )
    panaroo_160 = environment_prefix(
        "EC_K_PANAROO_160_RUN",
        "panaroo_1.6.0",
    )

    try:
        for path in official_gffs:
            shutil.copy2(path, inputs / path.name)

        private_name = f"G{args.group}_{locus}_{accession}.gff"
        private_input = inputs / private_name
        normalise_command: list[object] = [
            *builder_run,
            "python",
            root / "private_builder" / "normalise_private_gff.py",
            input_gff,
            private_input,
            "--group",
            args.group,
            "--locus",
            locus,
            "--accession",
            accession,
        ]
        if args.allow_atypical_boundaries:
            normalise_command.append("--allow-atypical-boundaries")
        run(normalise_command, work)
        run(
            ["bash", scripts / "edit_novel_gff.sh", private_input],
            work,
        )

        combined_gffs = sorted(inputs.glob("*.gff"))
        panaroo_output = work / "panaroo"
        run(
            [
                *panaroo_152,
                "panaroo",
                "--clean-mode",
                "sensitive",
                "--search_radius",
                "30000",
                "--trailing_recursive",
                "0",
                "--min_trailing_support",
                "1",
                "-t",
                args.threads,
                "-i",
                *combined_gffs,
                "-o",
                panaroo_output,
            ],
            work,
        )

        reconciliation_report = work / "reconciliation.tsv"
        reconciled_graph = work / "final_graph_reconciled.gml"
        run(
            [
                *builder_run,
                "python",
                root / "private_builder" / "reconcile_graph.py",
                reference_gffs,
                panaroo_output / "final_graph.gml",
                panaroo_output / "gene_data.csv",
                reconciled_graph,
                reconciliation_report,
                private_input.stem,
            ],
            work,
        )
        shutil.copy2(
            panaroo_output / "final_graph.gml",
            panaroo_output / "final_graph_before_reconciliation.gml",
        )
        shutil.copy2(reconciled_graph, panaroo_output / "final_graph.gml")

        generate_inputs = work / "gffs_for_generate"
        run(
            [
                *builder_run,
                "python",
                scripts / "prepare_gffs_for_generate.py",
                inputs,
                generate_inputs,
            ],
            work,
        )
        run(
            [
                *panaroo_160,
                "panaroo-generate-gffs",
                "-i",
                *sorted(generate_inputs.glob("*.gff")),
                "-o",
                panaroo_output,
                "-t",
                args.threads,
            ],
            work,
        )

        final_gffs = work / "final_gffs"
        run(
            [
                *builder_run,
                "python",
                scripts / "finalise_generated_gffs.py",
                panaroo_output / "postpanaroo_gffs",
                inputs,
                final_gffs,
                "--remove-private-short-hypotheticals",
            ],
            work,
        )

        for official_gff in official_gffs:
            shutil.copy2(
                official_gff,
                final_gffs / official_gff.name,
            )
        for official_gff in official_gffs:
            if sha256(official_gff) != sha256(
                final_gffs / official_gff.name
            ):
                raise RuntimeError(
                    f"Frozen reference GFF changed: {official_gff.name}"
                )
        if not (final_gffs / private_name).is_file():
            raise FileNotFoundError(final_gffs / private_name)

        metadata_path = work / "locus_metadata.csv"
        private_metadata(
            base_metadata,
            metadata_path,
            locus,
            phenotype,
            note,
        )

        database_name = f"EC-K-typing_private_{locus}.gbk"
        database_path = work / database_name
        run(
            [
                *builder_run,
                "python",
                scripts / "K-gff_to_gbk.py",
                *sorted(final_gffs.glob("*.gff")),
                "--metadata",
                metadata_path,
                "--output",
                database_path,
            ],
            work,
        )

        extracted = work / "references"
        run(
            [
                *builder_run,
                "kaptive",
                "extract",
                database_path,
                "--fna",
                extracted,
            ],
            work,
        )
        result_path = work / "self_test.tsv"
        run(
            [
                *builder_run,
                "kaptive",
                "assembly",
                database_path,
                *sorted(extracted.glob("*.fna")),
                "-o",
                result_path,
            ],
            work,
        )

        summary = validate_self_test(
            result_path,
            len(combined_gffs),
            locus,
            phenotype,
        )
        tool_versions = {
            "kaptive": capture(
                [*builder_run, "kaptive", "--version"],
                work,
            ),
            "panaroo_clustering": capture(
                [*panaroo_152, "panaroo", "--version"],
                work,
            ),
            "panaroo_gff_generation": capture(
                [*panaroo_160, "panaroo", "--version"],
                work,
            ),
            "panaroo_clustering_runtime": capture(
                [
                    *panaroo_152,
                    "python",
                    "-c",
                    (
                        "import Bio,numpy; "
                        "print('Biopython='+Bio.__version__, "
                        "'NumPy='+numpy.__version__)"
                    ),
                ],
                work,
            ),
        }

        final_database = output_dir / database_name
        final_results = output_dir / "validation.tsv"
        final_private_gff = output_dir / private_name
        final_report = output_dir / "reconciliation.tsv"
        shutil.copy2(database_path, final_database)
        shutil.copy2(result_path, final_results)
        shutil.copy2(final_gffs / private_name, final_private_gff)
        shutil.copy2(reconciliation_report, final_report)

        manifest = {
            "builder": "EC-K-Typing private database builder",
            "builder_version": BUILDER_VERSION,
            "created_utc": datetime.now(timezone.utc).isoformat(),
            "official_database_version": reference_version(config),
            "official_reference_records": len(official_gffs),
            "official_reference_gffs_preserved": True,
            "official_reference_sha256": reference_sha256(
                [
                    *official_gffs,
                    base_metadata,
                    config,
                    scripts / "K-gff_to_gbk.py",
                ],
                root,
            ),
            "private_locus": locus,
            "private_group": int(args.group),
            "private_accession": accession,
            "private_phenotype": phenotype,
            "atypical_boundaries_allowed": (
                args.allow_atypical_boundaries
            ),
            "input_gff_sha256": sha256(input_gff),
            "database_file": final_database.name,
            "database_sha256": sha256(final_database),
            "tool_versions": tool_versions,
            "validation": summary,
            "source_repository_modified": False,
        }
        (output_dir / "manifest.json").write_text(
            json.dumps(manifest, indent=2, sort_keys=True) + "\n",
            encoding="utf-8",
        )

        if not args.keep_work:
            shutil.rmtree(work)

        print("")
        print("Private database built successfully")
        print(f"Database: {final_database}")
        print(f"Validation: {final_results}")
        print(f"Manifest: {output_dir / 'manifest.json'}")
    except Exception as error:
        (output_dir / "FAILED.txt").write_text(
            f"{type(error).__name__}: {error}\n"
            f"Intermediates retained in: {work}\n",
            encoding="utf-8",
        )
        raise


if __name__ == "__main__":
    main()
