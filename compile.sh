#!/bin/bash
set -euo pipefail

GPU_INDEX=${GPU_INDEX:-0}
if [ -n "${GPU_ARCH:-}" ]; then
    CC="${GPU_ARCH}"
else
    if command -v nvidia-smi >/dev/null 2>&1; then
        CC=$(nvidia-smi --query-gpu=compute_cap --format=csv,noheader -i "${GPU_INDEX}" 2>/dev/null | head -n1 | tr -d '.')
    fi
    CC=${CC:-86} 
fi

VELOCITY_SET=${1:-}
ID=${2:-}

if [ -z "$VELOCITY_SET" ] || [ -z "$ID" ]; then
    echo "Usage: ./compile.sh <VELOCITY_SET> <ID>"
    echo "Example: ./compile.sh D3Q19 000"
    exit 1
fi

if [ "$VELOCITY_SET" != "D3Q19" ] && [ "$VELOCITY_SET" != "D3Q27" ]; then
    echo "Invalid VELOCITY_SET. Use 'D3Q19' or 'D3Q27'."
    exit 1
fi

SCRIPT_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" &>/dev/null && pwd)"

if [ -d "${SCRIPT_DIR}/src" ]; then
    BASE_DIR="${SCRIPT_DIR}"
elif [ -d "${SCRIPT_DIR}/../src" ]; then
    BASE_DIR="$(cd "${SCRIPT_DIR}/.." && pwd)"
else
    echo "Error: could not locate project root (missing 'src/' next to or above compile.sh)."
    echo "Checked: ${SCRIPT_DIR}/src and ${SCRIPT_DIR}/../src"
    exit 1
fi

SRC_DIR="${BASE_DIR}/src"
OUTPUT_DIR="${BASE_DIR}/bin/${VELOCITY_SET}"
EXECUTABLE="${OUTPUT_DIR}/${ID}sim_${VELOCITY_SET}_sm${CC}"

mkdir -p "${OUTPUT_DIR}"

echo "Project root detected: ${BASE_DIR}" 
echo "Compiling to ${EXECUTABLE}..."

nvcc -O3 --restrict \
     -gencode arch=compute_${CC},code=sm_${CC} \
     -gencode arch=compute_${CC},code=compute_${CC} \
     -rdc=true --ptxas-options=-v \
     --use_fast_math --fmad=true \
     -I"${SRC_DIR}" \
     -std=c++20 "${SRC_DIR}/main.cu" \
     -D${VELOCITY_SET} \
     -DENABLE_FP16=1 \
     --extended-lambda \
     -o "${EXECUTABLE}"

echo "Compilation completed successfully: ${EXECUTABLE}"
