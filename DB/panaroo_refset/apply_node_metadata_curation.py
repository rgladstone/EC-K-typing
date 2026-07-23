#!/usr/bin/env python3

"""Apply an explicit, complete node-metadata mapping to a Panaroo graph."""

import argparse
import csv
from collections import Counter
from copy import deepcopy
from pathlib import Path

import networkx as nx


EDITABLE_ATTRIBUTES = {"name", "annotation", "description"}
REQUIRED_COLUMNS = {
    "current_name",
    "desired_name",
    "desired_annotation",
    "desired_description",
    "reason",
    "historical_match_count",
}


def parse_arguments() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Replace the name, annotation and description of every Panaroo "
            "graph node using an explicitly curated CSV mapping."
        )
    )
    parser.add_argument("input_graph", type=Path)
    parser.add_argument("mapping_csv", type=Path)
    parser.add_argument("output_graph", type=Path)
    return parser.parse_args()


def load_mapping(path: Path) -> tuple[dict[str, dict[str, str]], Counter]:
    with path.open(newline="", encoding="utf-8-sig") as handle:
        reader = csv.DictReader(handle)
        fieldnames = set(reader.fieldnames or [])

        if fieldnames != REQUIRED_COLUMNS:
            raise ValueError(
                f"{path} must contain exactly these columns: "
                + ", ".join(sorted(REQUIRED_COLUMNS))
            )

        rows = list(reader)

    mapping: dict[str, dict[str, str]] = {}
    desired_names: set[str] = set()
    reasons: Counter = Counter()

    for row in rows:
        current_name = row["current_name"].strip()
        desired_name = row["desired_name"].strip()
        desired_annotation = row["desired_annotation"].strip()
        desired_description = row["desired_description"].strip()
        reason = row["reason"].strip()

        if not all(
            (
                current_name,
                desired_name,
                desired_annotation,
                desired_description,
                reason,
            )
        ):
            raise ValueError(f"Empty required mapping field in row: {row}")

        if current_name in mapping:
            raise ValueError(f"Duplicate current_name: {current_name}")

        if desired_name in desired_names:
            raise ValueError(f"Duplicate desired_name: {desired_name}")

        if "~~~" in desired_name:
            raise ValueError(
                f"Merged delimiter remains in desired_name: {desired_name}"
            )

        if ";" in desired_description:
            raise ValueError(
                f"Semicolon remains in desired_description for {current_name}"
            )

        desired_names.add(desired_name)
        reasons[reason] += 1
        mapping[current_name] = {
            "name": desired_name,
            "annotation": desired_annotation,
            "description": desired_description,
        }

    return mapping, reasons


def graph_names(graph: nx.Graph) -> dict[str, object]:
    names: dict[str, object] = {}

    for node, data in graph.nodes(data=True):
        name = str(data.get("name", ""))
        if not name:
            raise ValueError(f"Graph node has no name: {node}")
        if name in names:
            raise ValueError(f"Duplicate graph node name: {name}")
        names[name] = node

    return names


def without_editable_attributes(data: dict) -> dict:
    return {
        key: deepcopy(value)
        for key, value in data.items()
        if key not in EDITABLE_ATTRIBUTES
    }


def main() -> None:
    args = parse_arguments()

    if not args.input_graph.is_file():
        raise FileNotFoundError(args.input_graph)

    if not args.mapping_csv.is_file():
        raise FileNotFoundError(args.mapping_csv)

    if args.output_graph.exists():
        raise FileExistsError(
            f"Refusing to overwrite existing output: {args.output_graph}"
        )

    mapping, reasons = load_mapping(args.mapping_csv)

    # Use NetworkX's default label handling. Loading with label=None changes
    # the original GML label attributes when the graph is written again.
    graph = nx.read_gml(args.input_graph)
    names = graph_names(graph)

    missing = set(names) - set(mapping)
    extra = set(mapping) - set(names)

    if missing:
        raise ValueError(
            "Graph nodes missing from metadata mapping: "
            + ", ".join(sorted(missing))
        )

    if extra:
        raise ValueError(
            "Metadata mappings absent from graph: "
            + ", ".join(sorted(extra))
        )

    original_nodes = {
        node: without_editable_attributes(data)
        for node, data in graph.nodes(data=True)
    }
    original_edges = {
        (source, target): deepcopy(data)
        for source, target, data in graph.edges(data=True)
    }

    for current_name, node in names.items():
        curated = mapping[current_name]
        graph.nodes[node].update(curated)

    final_names = [
        str(data.get("name", ""))
        for _, data in graph.nodes(data=True)
    ]
    final_annotations = [
        str(data.get("annotation", ""))
        for _, data in graph.nodes(data=True)
    ]

    if len(final_names) != len(set(final_names)):
        raise RuntimeError("Curated graph node names are not unique")

    if any(name in {"", "No_name"} for name in final_names):
        raise RuntimeError("Curated graph contains an unnamed node")

    if any(name in {"", "No_name"} for name in final_annotations):
        raise RuntimeError("Curated graph contains an unnamed annotation")

    if any("~~~" in name for name in final_names):
        raise RuntimeError("Curated graph contains a merged cluster name")

    nx.write_gml(graph, args.output_graph)
    check = nx.read_gml(args.output_graph)

    if set(check.nodes) != set(graph.nodes):
        raise RuntimeError("Graph node identifiers changed during writing")

    if set(check.edges) != set(graph.edges):
        raise RuntimeError("Graph edges changed during writing")

    for node, data in check.nodes(data=True):
        if without_editable_attributes(data) != original_nodes[node]:
            raise RuntimeError(
                f"Unexpected non-metadata change in graph node {node}"
            )

    for source, target, data in check.edges(data=True):
        if data != original_edges[(source, target)]:
            raise RuntimeError(
                f"Unexpected edge change: {source} -- {target}"
            )

    print(
        f"Wrote {args.output_graph}: "
        f"{check.number_of_nodes()} nodes, "
        f"{check.number_of_edges()} edges"
    )
    print(f"Curated nodes: {len(mapping)}")

    for reason, count in sorted(reasons.items()):
        print(f"{reason}: {count}")


if __name__ == "__main__":
    main()
