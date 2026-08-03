#!/bin/bash
#SBATCH -p amd_m9_768
#SBATCH -N 1
#SBATCH -n 256
#SBATCH --mem=400G
#SBATCH -J euler_lg_pf_eta3_lr
#SBATCH -t 72:00:00
#SBATCH -o logs/sim_%j.log
#SBATCH -e logs/sim_%j.err

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
CASE_DIR="${CASE_DIR:-${SLURM_SUBMIT_DIR:-$SCRIPT_DIR}}"
CASE_NAME="${CASE_NAME:-pcon_false_eta3_local_refinement}"
cd "$CASE_DIR"
mkdir -p logs output

source /publicfs10/fs10-share/soft/share-soft/modules/module.sh
module purge
module load gcc/9.3.0
module load cmake/3.30.2

ROOT_DIR=${ROOT_DIR:-$HOME/SPHinXsys}
BUILD_DIR=${BUILD_DIR:-$ROOT_DIR/build_2d}

export Simbody_DIR="$ROOT_DIR/simbody/lib64/cmake/simbody"
export Eigen3_DIR="$ROOT_DIR/Eigen3/share/eigen3/cmake"
export TBB_DIR="$ROOT_DIR/yilai/TBB/lib64/cmake/TBB"
export spdlog_DIR="$ROOT_DIR/yilai/spdlog/lib64/cmake/spdlog"
export GTest_DIR="$ROOT_DIR/yilai/gtest/lib64/cmake/GTest"
export BOOST_ROOT="$ROOT_DIR/yilai/boost_1_82_0"

export LD_LIBRARY_PATH="$ROOT_DIR/simbody/lib64:$ROOT_DIR/yilai/TBB/lib:$ROOT_DIR/yilai/spdlog/lib:$ROOT_DIR/yilai/gtest/lib64:$ROOT_DIR/yilai/boost_1_82_0/lib:${LD_LIBRARY_PATH:-}"

export OMP_NUM_THREADS=$SLURM_NTASKS
export TBB_NUM_THREADS=$SLURM_NTASKS
export OMP_PROC_BIND=true
export OMP_PLACES=cores

TARGET_NAME="test_2d_eulerian_flow_around_cylinder_LG"
REMAP_TARGET_NAME="${TARGET_NAME}_remap"

echo "=== CMake Configuration (2D Eulerian Cylinder LG) ==="
cmake --fresh -S "$ROOT_DIR" -B "$BUILD_DIR" \
  -Wno-dev \
  -DCMAKE_BUILD_TYPE=Release \
  -DSPHINXSYS_BUILD_TESTS=ON \
  -DSPHINXSYS_2D=ON \
  -DSPHINXSYS_3D=OFF \
  -DSPHINXSYS_BUILD_2D_EXAMPLES=ON \
  -DSPHINXSYS_2D_EXAMPLE_ONLY="$TARGET_NAME" \
  -DSPHINXSYS_BUILD_3D_EXAMPLES=OFF \
  -DSPHINXSYS_BUILD_OPTIMIZATION_EXAMPLES=OFF \
  -DSPHINXSYS_BUILD_USER_EXAMPLES=OFF \
  -DSPHINXSYS_BUILD_MODULES=OFF \
  -DSPHINXSYS_BUILD_PYTHON_INTERFACE=OFF \
  -DSPHINXSYS_BUILD_UNIT_TESTS=OFF \
  -DSimbody_DIR="$Simbody_DIR" \
  -DEigen3_DIR="$Eigen3_DIR" \
  -DTBB_DIR="$TBB_DIR" \
  -Dspdlog_DIR="$spdlog_DIR" \
  -DGTest_DIR="$GTest_DIR" \
  -DBOOST_ROOT="$BOOST_ROOT"

echo "=== Building $TARGET_NAME and $REMAP_TARGET_NAME ==="
cmake --build "$BUILD_DIR" --target "$TARGET_NAME" "$REMAP_TARGET_NAME" --config Release --parallel 64

SHARED_BINARY="$BUILD_DIR/tests/2d_examples/test_2d_eulerian_flow_around_cylinder_LG/bin/test_2d_eulerian_flow_around_cylinder_LG"
SHARED_REMAP_BINARY="${SHARED_BINARY}_remap"
CASE_BINARY="$CASE_DIR/test_2d_eulerian_flow_around_cylinder_LG_${CASE_NAME}"
CASE_REMAP_BINARY="${CASE_BINARY}_remap"
cp -f "$SHARED_BINARY" "$CASE_BINARY"
cp -f "$SHARED_REMAP_BINARY" "$CASE_REMAP_BINARY"
chmod +x "$CASE_BINARY"
chmod +x "$CASE_REMAP_BINARY"

echo "=== Running simulation at $(date) ==="
srun -n 1 "$CASE_BINARY"

echo "=== Job completed at $(date) ==="
echo "=== Submitting post-processing job at $(date) ==="
sbatch ../../post.sh
