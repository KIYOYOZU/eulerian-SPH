#!/bin/bash
#SBATCH -p amd_m9_768
#SBATCH -N 1
#SBATCH -n 64
#SBATCH --mem=128G
#SBATCH -J eulerian_cylinder_2d
#SBATCH -t 24:00:00
#SBATCH -o logs/2d_cylinder_%j.log
#SBATCH -e logs/2d_cylinder_%j.err

set -euo pipefail

CASE_ROOT_DIR="${SLURM_SUBMIT_DIR:-/publicfs10/fs10-m9/home/m9s001530/eulerian-SPH/cases/2d_eulerian_flow_around_cylinder_LG}"
ROOT_DIR="${ROOT_DIR:-/publicfs10/fs10-m9/home/m9s001530/eulerian-SPH}"
SPHINXSYS_DEPS_ROOT="${SPHINXSYS_DEPS_ROOT:-/publicfs10/fs10-m9/home/m9s001530/SPHinXsys}"
BUILD_DIR="${BUILD_DIR:-$ROOT_DIR/build}"
CASE_TARGET_NAME="2d_eulerian_flow_around_cylinder_LG"
CASE_BIN_DIR="$BUILD_DIR/cases/$CASE_TARGET_NAME/bin"
RUN_ARGS="${RUN_ARGS:-}"

cd "$CASE_ROOT_DIR"
mkdir -p logs output

source /publicfs10/fs10-share/soft/share-soft/modules/module.sh
module purge >/dev/null 2>&1
module load gcc/9.3.0 >/dev/null 2>&1
module load cmake/3.30.2 >/dev/null 2>&1
echo "[ENV] Modules loaded."

export Eigen3_DIR="${Eigen3_DIR:-$SPHINXSYS_DEPS_ROOT/Eigen3/share/eigen3/cmake}"
export TBB_DIR="${TBB_DIR:-$SPHINXSYS_DEPS_ROOT/yilai/TBB/lib64/cmake/TBB}"
export BOOST_ROOT="${BOOST_ROOT:-$SPHINXSYS_DEPS_ROOT/yilai/boost_1_82_0}"
export BOOST_LIBRARYDIR="${BOOST_LIBRARYDIR:-$BOOST_ROOT/lib}"

export LD_LIBRARY_PATH="$SPHINXSYS_DEPS_ROOT/yilai/TBB/lib:$BOOST_LIBRARYDIR:${LD_LIBRARY_PATH:-}"

export OMP_NUM_THREADS="${OMP_NUM_THREADS:-$SLURM_NTASKS}"
export TBB_NUM_THREADS="${TBB_NUM_THREADS:-$SLURM_NTASKS}"
export OMP_PROC_BIND=true
export OMP_PLACES=cores

echo "=== CMake Configuration ($CASE_TARGET_NAME) ==="
cmake -S "$ROOT_DIR" -B "$BUILD_DIR" \
  -Wno-dev \
  -DCMAKE_BUILD_TYPE=Release \
  -DEigen3_DIR="$Eigen3_DIR" \
  -DTBB_DIR="$TBB_DIR" \
  -DBOOST_ROOT="$BOOST_ROOT" \
  -DBOOST_LIBRARYDIR="$BOOST_LIBRARYDIR"

echo "=== Building $CASE_TARGET_NAME ==="
cmake --build "$BUILD_DIR" --target "$CASE_TARGET_NAME" --config Release -j 64

echo "=== Running $CASE_TARGET_NAME at $(date) ==="
if [ -n "$RUN_ARGS" ]; then
  echo "[RUN] srun -n 1 \"$CASE_BIN_DIR/$CASE_TARGET_NAME\" $RUN_ARGS"
  srun -n 1 "$CASE_BIN_DIR/$CASE_TARGET_NAME" $RUN_ARGS
else
  echo "[RUN] srun -n 1 \"$CASE_BIN_DIR/$CASE_TARGET_NAME\""
  srun -n 1 "$CASE_BIN_DIR/$CASE_TARGET_NAME"
fi

echo "=== Job completed at $(date) ==="
