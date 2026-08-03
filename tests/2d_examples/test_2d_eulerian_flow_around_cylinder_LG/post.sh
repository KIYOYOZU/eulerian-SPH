#!/bin/bash
#SBATCH -p amd_m9_768
#SBATCH -N 1
#SBATCH -n 64
#SBATCH --mem=100G
#SBATCH -J euler_lg_post
#SBATCH -t 02:00:00
#SBATCH -o logs/post_%j.log
#SBATCH -e logs/post_%j.err

CASE_DIR="${CASE_DIR:-${SLURM_SUBMIT_DIR:-$(pwd)}}"

if [[ -f "$CASE_DIR/postprocess_cl_cd.py" ]]; then
    SHARED_DIR="$(cd "$CASE_DIR" && pwd)"
elif [[ -f "$CASE_DIR/../postprocess_cl_cd.py" ]]; then
    SHARED_DIR="$(cd "$CASE_DIR/.." && pwd)"
elif [[ -f "$CASE_DIR/../../postprocess_cl_cd.py" ]]; then
    SHARED_DIR="$(cd "$CASE_DIR/../.." && pwd)"
else
    echo "[ERROR] Cannot locate postprocess_cl_cd.py from CASE_DIR=$CASE_DIR"
    exit 1
fi

cd "$CASE_DIR"
mkdir -p logs results

source /publicfs10/fs10-share/soft/share-soft/modules/module.sh
module purge
module load miniforge/25.9 || {
    echo "[ERROR] miniforge module not found"
    exit 1
}

python3 -c "import numpy, matplotlib" 2>/dev/null || \
    python3 -m pip install --user -q -i https://pypi.tuna.tsinghua.edu.cn/simple numpy matplotlib

echo "=== Post-processing started at $(date) ==="
echo "CASE_DIR=$CASE_DIR  SHARED_DIR=$SHARED_DIR"

srun -n 1 python3 -u "$SHARED_DIR/postprocess_cl_cd.py" \
    --case-dir "$CASE_DIR" \
    --output-dir output \
    --results-dir results

echo "=== Post-processing completed at $(date) ==="
echo ""
echo "Generated files:"
ls -lh results/ 2>/dev/null || echo "  No files in results/"
