#!/usr/bin/env python3

"""Restore frozen EC-K-Typing metadata in a Panaroo graph."""

import argparse
import csv
import hashlib
from collections import Counter
from copy import deepcopy
from pathlib import Path

import networkx as nx


EDITABLE_ATTRIBUTES = {"name", "annotation", "description"}


def parse_arguments() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Reconcile a Panaroo graph against frozen reference GFF feature "
            "metadata while assigning local names to merge/split artefacts."
        )
    )
    parser.add_argument("reference_gff_dir", type=Path)
    parser.add_argument("input_graph", type=Path)
    parser.add_argument("gene_data", type=Path)
    parser.add_argument("output_graph", type=Path)
    parser.add_argument("report", type=Path)
    parser.add_argument("private_gff_stem")
    return parser.parse_args()


def parse_attributes(value: str) -> dict[str, str]:
    return {
        key: item
        for field in value.split(";")
        if "=" in field
        for key, item in [field.split("=", 1)]
    }


def load_reference_metadata(
    directory: Path,
) -> dict[tuple[str, str], tuple[str, str, str]]:
    if not directory.is_dir():
        raise NotADirectoryError(directory)

    metadata: dict[tuple[str, str], tuple[str, str, str]] = {}
    gffs = sorted(directory.glob("*.gff"))
    if not gffs:
        raise ValueError(f"No reference GFFs found in {directory}")

    for path in gffs:
        stem = path.stem
        for line_number, line in enumerate(
            path.read_text(encoding="utf-8").splitlines(),
            1,
        ):
            if line == "##FASTA":
                break
            if not line or line.startswith("#"):
                continue
            fields = line.split("\t")
            if len(fields) != 9 or fields[2] not in {
                "CDS",
                "candidate_gene",
            }:
                continue

            attributes = parse_attributes(fields[8])
            annotation_id = attributes.get(
                "ID", attributes.get("locus_tag", "")
            )
            name = attributes.get("name", attributes.get("Name", ""))
            description = attributes.get(
                "description", attributes.get("product", "")
            )

            if not annotation_id or not name:
                raise ValueError(
                    f"Missing ID or name in {path}:{line_number}"
                )

            key = (stem, annotation_id)
            value = (name, name, description)
            if key in metadata and metadata[key] != value:
                raise ValueError(
                    f"Conflicting metadata for {key}: "
                    f"{metadata[key]} versus {value}"
                )
            metadata[key] = value

    return metadata


def load_gene_data(path: Path) -> dict[str, dict[str, str]]:
    if not path.is_file():
        raise FileNotFoundError(path)

    with path.open(newline="", encoding="utf-8") as handle:
        rows = list(csv.DictReader(handle))

    required = {
        "gff_file",
        "clustering_id",
        "annotation_id",
        "prot_sequence",
        "gene_name",
        "description",
    }
    if not rows or not required.issubset(rows[0]):
        raise ValueError(
            f"{path} is missing required columns: "
            + ", ".join(sorted(required))
        )

    result: dict[str, dict[str, str]] = {}
    for row in rows:
        clustering_id = row["clustering_id"]
        if clustering_id in result:
            raise ValueError(
                f"Duplicate clustering_id in {path}: {clustering_id}"
            )
        result[clustering_id] = row
    return result


def node_sequence_ids(data: dict[str, object]) -> list[str]:
    value = data.get("seqIDs", [])
    if isinstance(value, str):
        return [value]
    return [str(item) for item in value]


def private_metadata(
    rows: list[dict[str, str]],
) -> tuple[str, str, str]:
    proteins = sorted(
        {
            row.get("prot_sequence", "")
            for row in rows
            if row.get("prot_sequence", "")
        }
    )
    seed = "|".join(proteins)
    if not seed:
        dna_sequences = sorted(
            {
                row.get("dna_sequence", "")
                for row in rows
                if row.get("dna_sequence", "")
            }
        )
        seed = "|".join(dna_sequences)
    if not seed:
        raise ValueError(
            "Cannot assign a deterministic name to a private cluster "
            "without a protein or DNA sequence"
        )
    digest = hashlib.sha256(seed.encode("utf-8")).hexdigest()[:12]
    name = f"private_group_{digest}"

    descriptions = [
        row.get("description", "").strip()
        for row in rows
        if row.get("description", "").strip()
    ]
    description = (
        Counter(descriptions).most_common(1)[0][0]
        if descriptions
        else "hypothetical protein"
    )
    description = " ".join(description.replace(";", ",").split())
    return name, name, description


def composite_metadata(
    matches: set[tuple[str, str, str]],
) -> tuple[str, str, str]:
    names = sorted(item[0] for item in matches)
    digest = hashlib.sha256(
        "|".join(names).encode("utf-8")
    ).hexdigest()[:12]
    name = f"private_composite_{digest}"
    description = "Local composite of historical clusters: " + ", ".join(
        names
    )
    return name, name, description


def split_metadata(
    historical_name: str,
    rows: list[dict[str, str]],
) -> tuple[str, str, str]:
    proteins = sorted(
        row.get("prot_sequence", "")
        for row in rows
        if row.get("prot_sequence", "")
    )
    dna_sequences = sorted(
        row.get("dna_sequence", "")
        for row in rows
        if row.get("dna_sequence", "")
    )
    clustering_ids = sorted(
        row.get("clustering_id", "")
        for row in rows
        if row.get("clustering_id", "")
    )
    seed = "|".join(
        [
            historical_name,
            *proteins,
            *dna_sequences,
            *clustering_ids,
        ]
    )
    digest = hashlib.sha256(seed.encode("utf-8")).hexdigest()[:12]
    name = f"private_split_{digest}"

    descriptions = [
        row.get("description", "").strip()
        for row in rows
        if row.get("description", "").strip()
    ]
    description = (
        Counter(descriptions).most_common(1)[0][0]
        if descriptions
        else f"Local split of historical cluster {historical_name}"
    )
    description = " ".join(description.replace(";", ",").split())
    return name, name, description


def main() -> None:
    args = parse_arguments()

    if args.output_graph.exists():
        raise FileExistsError(args.output_graph)
    if args.report.exists():
        raise FileExistsError(args.report)

    reference = load_reference_metadata(args.reference_gff_dir)
    gene_data = load_gene_data(args.gene_data)
    # Use NetworkX's default label handling so the original GML label values
    # remain the node identifiers when the graph is written again.
    graph = nx.read_gml(args.input_graph)
    original = graph.copy()

    assignments: dict[object, tuple[str, str, str]] = {}
    contexts: dict[object, dict[str, object]] = {}
    problems: list[str] = []

    for node, data in graph.nodes(data=True):
        sequence_ids = node_sequence_ids(data)
        node_rows = [
            gene_data[sequence_id]
            for sequence_id in sequence_ids
            if sequence_id in gene_data
        ]
        if not node_rows:
            problems.append(f"Node {node} has no gene_data rows")
            continue

        matches = {
            reference[(row["gff_file"], row["annotation_id"])]
            for row in node_rows
            if (row["gff_file"], row["annotation_id"]) in reference
        }

        if len(matches) > 1:
            desired = composite_metadata(matches)
            status = "local_historical_composite"
        elif len(matches) == 1:
            desired = next(iter(matches))
            status = "restore_historical"
        else:
            desired = private_metadata(node_rows)
            status = "new_private_cluster"

        assignments[node] = desired
        contexts[node] = {
            "status": status,
            "current_name": str(data.get("name", "")),
            "rows": node_rows,
            "members": len(sequence_ids),
            "private_members": sum(
                row["gff_file"] == args.private_gff_stem
                for row in node_rows
            ),
            "historical_names": ",".join(
                sorted(item[0] for item in matches)
            ),
        }

    desired_nodes: dict[str, list[object]] = {}
    for node, desired in assignments.items():
        desired_nodes.setdefault(desired[0], []).append(node)
    for name, nodes in desired_nodes.items():
        if not name or len(nodes) == 1:
            continue

        winner = max(
            nodes,
            key=lambda node: (
                int(contexts[node]["private_members"]),
                int(contexts[node]["members"]),
                str(node),
            ),
        )
        for node in nodes:
            if node == winner:
                continue
            assignments[node] = split_metadata(
                name,
                contexts[node]["rows"],
            )
            contexts[node]["status"] = (
                f"{contexts[node]['status']}+local_historical_split"
            )

    final_names = [desired[0] for desired in assignments.values()]
    if len(final_names) != len(set(final_names)):
        problems.append(
            "Local reconciliation produced duplicate graph names"
        )

    rows_for_report = [
        {
            "node": str(node),
            "status": str(contexts[node]["status"]),
            "current_name": str(contexts[node]["current_name"]),
            "desired_name": assignments[node][0],
            "description": assignments[node][2],
            "members": str(contexts[node]["members"]),
            "private_members": str(contexts[node]["private_members"]),
            "historical_names": str(contexts[node]["historical_names"]),
        }
        for node in graph.nodes
        if node in assignments
    ]

    args.report.parent.mkdir(parents=True, exist_ok=True)
    with args.report.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=[
                "node",
                "status",
                "current_name",
                "desired_name",
                "description",
                "members",
                "private_members",
                "historical_names",
            ],
            delimiter="\t",
        )
        writer.writeheader()
        writer.writerows(rows_for_report)

    if problems:
        raise RuntimeError(
            "Graph reconciliation failed:\n"
            + "\n".join(f"- {problem}" for problem in problems)
            + f"\nReport: {args.report}"
        )

    for node, desired in assignments.items():
        name, annotation, description = desired
        graph.nodes[node]["name"] = name
        graph.nodes[node]["annotation"] = annotation
        graph.nodes[node]["description"] = description

    if set(original.nodes) != set(graph.nodes):
        raise RuntimeError("Graph nodes changed during reconciliation")
    if set(original.edges) != set(graph.edges):
        raise RuntimeError("Graph edges changed during reconciliation")

    for node in original.nodes:
        before = {
            key: value
            for key, value in original.nodes[node].items()
            if key not in EDITABLE_ATTRIBUTES
        }
        after = {
            key: value
            for key, value in graph.nodes[node].items()
            if key not in EDITABLE_ATTRIBUTES
        }
        if before != after:
            raise RuntimeError(
                f"Unexpected non-metadata change in node {node}"
            )

    args.output_graph.parent.mkdir(parents=True, exist_ok=True)
    nx.write_gml(graph, args.output_graph)
    check = nx.read_gml(args.output_graph)

    if set(check.nodes) != set(graph.nodes):
        raise RuntimeError(
            "Graph node identifiers changed during GML writing"
        )
    if set(check.edges) != set(graph.edges):
        raise RuntimeError("Graph edges changed during GML writing")

    for node, data in check.nodes(data=True):
        before = {
            key: deepcopy(value)
            for key, value in graph.nodes[node].items()
            if key not in EDITABLE_ATTRIBUTES
        }
        after = {
            key: deepcopy(value)
            for key, value in data.items()
            if key not in EDITABLE_ATTRIBUTES
        }
        if before != after:
            raise RuntimeError(
                f"Unexpected non-metadata change after writing node {node}"
            )

    for source, target, data in check.edges(data=True):
        if data != graph.edges[source, target]:
            raise RuntimeError(
                f"Unexpected edge change after writing: "
                f"{source} -- {target}"
            )

    new_count = sum(
        row["status"] == "new_private_cluster"
        for row in rows_for_report
    )
    composite_count = sum(
        "composite" in row["status"] for row in rows_for_report
    )
    split_count = sum(
        "split" in row["status"] for row in rows_for_report
    )
    print(f"Nodes: {check.number_of_nodes()}")
    print(f"Edges: {check.number_of_edges()}")
    print(f"New private clusters: {new_count}")
    print(f"Local composite clusters: {composite_count}")
    print(f"Local split clusters: {split_count}")
    print(f"Wrote graph: {args.output_graph}")
    print(f"Wrote report: {args.report}")


if __name__ == "__main__":
    main()
