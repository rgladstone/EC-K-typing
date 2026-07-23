#!/usr/bin/env python3

import argparse
import csv
from pathlib import Path

import networkx as nx


def parse_arguments():
    parser = argparse.ArgumentParser(
        description=(
            "Replace merged Panaroo cluster names and annotations using "
            "an explicitly curated CSV mapping."
        )
    )
    parser.add_argument("input_graph", type=Path)
    parser.add_argument("mapping_csv", type=Path)
    parser.add_argument("output_graph", type=Path)
    return parser.parse_args()


def load_mapping(path):
    required = {
        "merged_name",
        "curated_name",
        "curated_description",
    }

    with path.open(newline="") as handle:
        reader = csv.DictReader(handle)

        if set(reader.fieldnames or []) != required:
            raise ValueError(
                f"{path} must contain exactly these columns: "
                + ", ".join(sorted(required))
            )

        rows = list(reader)

    mapping = {}

    for row in rows:
        merged_name = row["merged_name"].strip()
        curated_name = row["curated_name"].strip()
        curated_description = row["curated_description"].strip()

        if not merged_name or not curated_name or not curated_description:
            raise ValueError(f"Empty mapping field in row: {row}")

        if merged_name in mapping:
            raise ValueError(f"Duplicate merged name: {merged_name}")

        mapping[merged_name] = {
            "name": curated_name,
            "annotation": curated_name,
            "description": curated_description,
        }

    return mapping


def main():
    args = parse_arguments()

    if not args.input_graph.is_file():
        raise FileNotFoundError(args.input_graph)

    if not args.mapping_csv.is_file():
        raise FileNotFoundError(args.mapping_csv)

    if args.output_graph.exists():
        raise FileExistsError(
            f"Refusing to overwrite existing output: {args.output_graph}"
        )

    mapping = load_mapping(args.mapping_csv)
    graph = nx.read_gml(args.input_graph)

    merged_nodes = {}

    for node, data in graph.nodes(data=True):
        name = str(data.get("name", ""))

        if "~~~" not in name:
            continue

        if name in merged_nodes:
            raise ValueError(f"Duplicate merged graph name: {name}")

        merged_nodes[name] = node

    missing = set(merged_nodes) - set(mapping)
    extra = set(mapping) - set(merged_nodes)

    if missing:
        raise ValueError(
            "Merged graph nodes missing from mapping: "
            + ", ".join(sorted(missing))
        )

    if extra:
        raise ValueError(
            "Mapping entries absent from graph: "
            + ", ".join(sorted(extra))
        )

    unchanged_names = {
        str(data.get("name", ""))
        for node, data in graph.nodes(data=True)
        if node not in set(merged_nodes.values())
    }

    collisions = {
        values["name"]
        for values in mapping.values()
        if values["name"] in unchanged_names
    }

    if collisions:
        raise ValueError(
            "Curated names collide with unchanged graph nodes: "
            + ", ".join(sorted(collisions))
        )

    node_count = graph.number_of_nodes()
    edge_count = graph.number_of_edges()

    for merged_name, node in merged_nodes.items():
        curated = mapping[merged_name]
        data = graph.nodes[node]

        data["name"] = curated["name"]
        data["annotation"] = curated["annotation"]
        data["description"] = curated["description"]

        print(f"{merged_name} -> {curated['name']}")

    remaining = [
        str(data.get("name", ""))
        for _, data in graph.nodes(data=True)
        if "~~~" in str(data.get("name", ""))
    ]

    if remaining:
        raise RuntimeError(
            "Merged names remain after curation: "
            + ", ".join(sorted(remaining))
        )

    nx.write_gml(graph, args.output_graph)

    check = nx.read_gml(args.output_graph)

    if check.number_of_nodes() != node_count:
        raise RuntimeError("Node count changed during GML writing")

    if check.number_of_edges() != edge_count:
        raise RuntimeError("Edge count changed during GML writing")

    remaining_after_write = [
        str(data.get("name", ""))
        for _, data in check.nodes(data=True)
        if "~~~" in str(data.get("name", ""))
    ]

    if remaining_after_write:
        raise RuntimeError(
            "Merged names remain in written graph: "
            + ", ".join(sorted(remaining_after_write))
        )

    print(
        f"Wrote {args.output_graph}: "
        f"{check.number_of_nodes()} nodes, "
        f"{check.number_of_edges()} edges, "
        f"{len(mapping)} curated clusters"
    )


if __name__ == "__main__":
    main()
