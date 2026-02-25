#!/usr/bin/env bash
set -euo pipefail

# Create/update BEES conda environment 
SCRIPT_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" &>/dev/null && pwd)"
REPO_ROOT="$(cd "$SCRIPT_DIR/.." && pwd)"
ENV_NAME=bees_env
ENV_FILE="$REPO_ROOT/environment.yaml"

if command -v micromamba &>/dev/null; then
  COMMAND_PKG=micromamba
elif command -v mamba &>/dev/null; then
  COMMAND_PKG=mamba
elif command -v conda &>/dev/null; then
  COMMAND_PKG=conda
else
  echo "Conda, Mamba, or Micromamba is required."
  exit 1
fi

echo "Using $COMMAND_PKG"

# Use CONDA_ALWAYS_YES to avoid prompts (some mamba versions reject -y/--yes for env commands)
export CONDA_ALWAYS_YES=true

if $COMMAND_PKG env list | awk '{print $1}' | sed 's/^\*//' | grep -Fxq "$ENV_NAME"; then
  echo "Updating environment: $ENV_NAME"
  $COMMAND_PKG env update -n "$ENV_NAME" -f "$ENV_FILE"
else
  echo "Creating environment: $ENV_NAME"
  $COMMAND_PKG env create -n "$ENV_NAME" -f "$ENV_FILE"
fi

echo "BEES environment ready. Activate with: conda activate $ENV_NAME"
