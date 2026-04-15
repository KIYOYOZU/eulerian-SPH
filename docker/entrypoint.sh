#!/usr/bin/env bash

set -euo pipefail

ROOT_DIR="${ROOT_DIR:-/opt/eulerian-SPH}"
BUILD_DIR="${BUILD_DIR:-$ROOT_DIR/build-docker}"
BUILD_TYPE="${BUILD_TYPE:-Release}"
CMAKE_GENERATOR="${CMAKE_GENERATOR:-Ninja}"
CASE_NAME="${CASE_NAME:-2d_eulerian_taylor_green_LG}"
BUILD_JOBS="${BUILD_JOBS:-}"
CONFIGURE_ONLY="${CONFIGURE_ONLY:-0}"
SKIP_BUILD="${SKIP_BUILD:-0}"
RUN_ARGS="${RUN_ARGS:-}"

declare -a SUPPORTED_CASES=(
  "2d_eulerian_taylor_green_LG"
  "2d_eulerian_flow_around_cylinder_LG"
  "2d_eulerian_supersonic_flow_around_cylinder"
)

case_supported=0
for candidate in "${SUPPORTED_CASES[@]}"; do
  if [[ "$candidate" == "$CASE_NAME" ]]; then
    case_supported=1
    break
  fi
done

if [[ "$case_supported" -ne 1 ]]; then
  echo "[ERROR] Unsupported CASE_NAME: $CASE_NAME" >&2
  echo "[ERROR] Supported cases: ${SUPPORTED_CASES[*]}" >&2
  exit 1
fi

CASE_DIR="$ROOT_DIR/cases/$CASE_NAME"
BINARY_PATH="$BUILD_DIR/cases/$CASE_NAME/bin/$CASE_NAME"

if [[ ! -d "$ROOT_DIR" ]]; then
  echo "[ERROR] ROOT_DIR does not exist: $ROOT_DIR" >&2
  exit 1
fi

if [[ ! -d "$CASE_DIR" ]]; then
  echo "[ERROR] Case directory does not exist: $CASE_DIR" >&2
  exit 1
fi

mkdir -p "$BUILD_DIR"

if [[ -z "$BUILD_JOBS" || ! "$BUILD_JOBS" =~ ^[0-9]+$ || "$BUILD_JOBS" -le 0 ]]; then
  BUILD_JOBS="$(nproc)"
fi

configure_project() {
  echo "[INFO] Configuring project"
  cmake -S "$ROOT_DIR" -B "$BUILD_DIR" \
    -G "$CMAKE_GENERATOR" \
    -DCMAKE_BUILD_TYPE="$BUILD_TYPE"
}

if [[ ! -f "$BUILD_DIR/CMakeCache.txt" ]]; then
  configure_project
elif [[ "${FORCE_CONFIGURE:-0}" == "1" ]]; then
  configure_project
fi

if [[ "$CONFIGURE_ONLY" == "1" ]]; then
  echo "[INFO] Configuration finished. CONFIGURE_ONLY=1, exiting."
  exit 0
fi

if [[ "$SKIP_BUILD" != "1" ]]; then
  echo "[INFO] Building target: $CASE_NAME"
  cmake --build "$BUILD_DIR" --target "$CASE_NAME" --parallel "$BUILD_JOBS"
fi

if [[ ! -x "$BINARY_PATH" ]]; then
  echo "[ERROR] Executable not found: $BINARY_PATH" >&2
  exit 1
fi

mkdir -p "$CASE_DIR/output" "$CASE_DIR/restart" "$CASE_DIR/reload" "$CASE_DIR/logs"
cd "$CASE_DIR"

declare -a runtime_args=()
if [[ $# -gt 0 ]]; then
  runtime_args=("$@")
elif [[ -n "$RUN_ARGS" ]]; then
  read -r -a runtime_args <<< "$RUN_ARGS"
fi

echo "[INFO] Running case: $CASE_NAME"
echo "[INFO] Working directory: $CASE_DIR"
echo "[INFO] Executable: $BINARY_PATH"

exec "$BINARY_PATH" "${runtime_args[@]}"
