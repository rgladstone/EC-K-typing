#!/bin/bash
#SBATCH --job-name=K_panaroo
#SBATCH --cpus-per-task=36
#SBATCH --mem=32G
#SBATCH --time=04:00:00
#SBATCH --output=DB/panaroo_refset/slurm_logs/panaroo_%j.out
#SBATCH --error=DB/panaroo_refset/slurm_logs/panaroo_%j.err

set -euo pipefail

cd "$SLURM_SUBMIT_DIR"

if ! command -v conda >/dev/null 2>&1; then
    echo "ERROR: conda is not available in the submitted job PATH" >&2
    exit 1
fi
conda_base=$(conda info --base)
source "$conda_base/etc/profile.d/conda.sh"
conda activate "${PANAROO_152_ENV:-panaroo_1.5.2}"
python -c '
from importlib.metadata import version
assert version("panaroo") == "1.5.2", version("panaroo")
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

echo "Panaroo executable: $(command -v panaroo)"
panaroo --version
echo "Input GFF count: $(find DB/panaroo_refset/gffs -maxdepth 1 -type f -name '*.gff' | wc -l)"

if [[ -e panaroo_add_novel ]]; then
    echo "ERROR: output path already exists: panaroo_add_novel" >&2
    exit 1
fi

panaroo \
    --clean-mode sensitive \
    --search_radius 30000 \
    --trailing_recursive 0 \
    --min_trailing_support 1 \
    -t "$SLURM_CPUS_PER_TASK" \
    -i DB/panaroo_refset/gffs/*.gff \
    -o panaroo_add_novel
