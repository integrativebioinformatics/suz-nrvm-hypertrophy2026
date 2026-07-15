#!/usr/bin/env bash
# Helper script to activate conda environments non-interactively
# Source this in any script that needs a specific conda environment
# Usage: source scripts/bash/_activate_env.sh env_qc_align
#
# This centralizes conda activation logic so it's not copy-pasted in every numbered script.

if [ -z "$1" ]; then
  echo "Error: _activate_env.sh requires an environment name as argument" >&2
  exit 1
fi

ENV_NAME="$1"

# Initialize conda shell hooks non-interactively
CONDA_BASE="$(conda info --base 2>/dev/null)" || {
  echo "Error: conda is not installed or not in PATH" >&2
  exit 1
}

# Source the conda shell hook (works in bash/sh/zsh/etc.)
if [ -f "$CONDA_BASE/etc/profile.d/conda.sh" ]; then
  source "$CONDA_BASE/etc/profile.d/conda.sh"
else
  echo "Error: conda shell hook not found at $CONDA_BASE/etc/profile.d/conda.sh" >&2
  exit 1
fi

# Activate the requested environment
conda activate "$ENV_NAME" 2>/dev/null || {
  echo "Error: Failed to activate conda environment '$ENV_NAME'" >&2
  echo "       Make sure the environment exists. Create it with:" >&2
  echo "       conda env create -f envs/environment_${ENV_NAME#env_}.yml" >&2
  exit 1
}

echo "[INFO] Activated conda environment: $ENV_NAME" >&2
