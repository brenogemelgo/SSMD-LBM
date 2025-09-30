#!/bin/bash
set -euo pipefail

if [ "$#" -ne 2 ]; then
    echo "Error: Usage: ./post.sh <velocity_set> <id>"
    echo "Example: ./post.sh D3Q19 000"
    exit 1
fi

VELOCITY_SET=$1
SIM_ID=$2

if [[ "$OSTYPE" == "linux-gnu"* ]]; then
    PYTHON_CMD="python3"
elif [[ "$OSTYPE" == "msys" || "$OSTYPE" == "cygwin" || "$OSTYPE" == "win32" ]]; then
    PYTHON_CMD="python"
else
    echo "Operating system not recognized. Trying python3 by default."
    PYTHON_CMD="python3"
fi

$PYTHON_CMD processSteps.py "$VELOCITY_SET" "$SIM_ID"
