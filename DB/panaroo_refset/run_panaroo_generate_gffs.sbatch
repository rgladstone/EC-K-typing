#!/bin/bash
#SBATCH --job-name=K_generate_gffs
#SBATCH --cpus-per-task=36
#SBATCH --mem=16G
#SBATCH --time=01:00:00
#SBATCH --output=DB/panaroo_refset/slurm_logs/generate_gffs_%j.out
#SBATCH --error=DB/panaroo_refset/slurm_logs/generate_gffs_%j.err

set -euo pipefail

cd "$SLURM_SUBMIT_DIR"

if ! command -v conda >/dev/null 2>&1; then
    echo "ERROR: conda is not available in the submitted job PATH" >&2
    exit 1
fi
conda_base=$(conda info --base)
source "$conda_base/etc/profile.d/conda.sh"
conda activate "${PANAROO_160_ENV:-panaroo_1.6.0}"

python -c '
from importlib.metadata import version
assert version("panaroo") == "1.6.0", version("panaroo")
print("Panaroo:", version("panaroo"))
'
python -c '
import Bio
import numpy
assert Bio.__version__ == "1.84", Bio.__version__
assert numpy.__version__ == "1.26.4", numpy.__version__
print("Biopython:", Bio.__version__)
print("NumPy:", numpy.__version__)
'

echo "Panaroo executable: $(command -v panaroo-generate-gffs)"
echo "Input GFF count: $(find panaroo_add_novel/gffs_for_generate -maxdepth 1 -type f -name '*.gff' | wc -l)"

source_count=$(
    find DB/panaroo_refset/gffs \
        -maxdepth 1 \
        -type f \
        -name '*.gff' |
    wc -l
)
generate_count=$(
    find panaroo_add_novel/gffs_for_generate \
        -maxdepth 1 \
        -type f \
        -name '*.gff' |
    wc -l
)

if [[ "$source_count" -ne "$generate_count" ]]; then
    echo "ERROR: source and generate-input GFF counts differ" >&2
    echo "source=$source_count generate=$generate_count" >&2
    exit 1
fi

if grep -Rql \
    --include='*.gff' \
    'panaroo_ID=' \
    panaroo_add_novel/gffs_for_generate; then
    echo "ERROR: legacy panaroo_ID fields remain in generate inputs" >&2
    exit 1
fi

if [[ -e panaroo_add_novel/postpanaroo_gffs ]]; then
    echo "ERROR: generated GFF output already exists" >&2
    echo "Move panaroo_add_novel/postpanaroo_gffs aside first." >&2
    exit 1
fi

python -c '
import networkx as nx
graph = nx.read_gml("panaroo_add_novel/final_graph.gml", label=None)
merged = [
    data.get("name", "")
    for _, data in graph.nodes(data=True)
    if "~~~" in str(data.get("name", ""))
]
assert not merged, merged
print("Curated graph:", graph.number_of_nodes(), "nodes,", graph.number_of_edges(), "edges")
'

panaroo-generate-gffs \
    -i panaroo_add_novel/gffs_for_generate/*.gff \
    -o panaroo_add_novel \
    -t "$SLURM_CPUS_PER_TASK"
