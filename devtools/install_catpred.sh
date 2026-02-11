#!/usr/bin/env bash
set -euo pipefail

# Clone CatPred, download pretrained data, create catpred env, write .env.bees.


SCRIPT_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" &>/dev/null && pwd)"
REPO_ROOT="$(cd "$SCRIPT_DIR/.." && pwd)"
CLONE_ROOT="$(dirname "$REPO_ROOT")"
CATPRED_REPO="${CLONE_ROOT}/CatPred"
CATPRED_PIPELINE="${CLONE_ROOT}/catpred_pipeline"
CATPRED_DATA_URL="${CATPRED_DATA_URL:-https://catpred.s3.us-east-1.amazonaws.com/capsule_data_update.tar.gz}"
CATPRED_GIT_URL="${CATPRED_FORK_URL:-https://github.com/DanaResearchGroup/CatPred.git}"

echo "CatPred clone: $CATPRED_GIT_URL"
echo "Clone root:    $CLONE_ROOT"
echo "CatPred path:  $CATPRED_REPO"

# Use same package manager as install_bees.sh (micromamba / mamba / conda)
if command -v micromamba &>/dev/null; then
  COMMAND_PKG=micromamba
elif command -v mamba &>/dev/null; then
  COMMAND_PKG=mamba
elif command -v conda &>/dev/null; then
  COMMAND_PKG=conda
else
  echo "Conda, Mamba, or Micromamba is required for CatPred env."
  exit 1
fi
echo "Using $COMMAND_PKG for CatPred env"

# Clone
if [[ ! -d "${CATPRED_REPO}/.git" ]]; then
  echo "Cloning CatPred into ${CATPRED_REPO}..."
  git clone "${CATPRED_GIT_URL}" "${CATPRED_REPO}"
else
  echo "CatPred already cloned at ${CATPRED_REPO}."
fi

# Pretrained data
mkdir -p "${CATPRED_PIPELINE}"
if [[ ! -f "${CATPRED_PIPELINE}/capsule_data_update.tar.gz" ]]; then
  echo "Downloading pretrained data..."
  (cd "${CATPRED_PIPELINE}" && (wget -q --show-progress "${CATPRED_DATA_URL}" || curl -L -o capsule_data_update.tar.gz "${CATPRED_DATA_URL}")) || true
fi
if [[ ! -d "${CATPRED_PIPELINE}/kcat" ]] && [[ ! -d "${CATPRED_PIPELINE}/production" ]]; then
  echo "Extracting pretrained data..."
  (cd "${CATPRED_PIPELINE}" && tar -xzf capsule_data_update.tar.gz 2>/dev/null) || true
fi

# Checkpoint base
CHECKPOINT_BASE="${CATPRED_PIPELINE}"
if [[ -d "${CATPRED_PIPELINE}/production" ]]; then
  CHECKPOINT_BASE="${CATPRED_PIPELINE}/production"
elif [[ -d "${CATPRED_PIPELINE}/capsule_data_update" ]]; then
  CHECKPOINT_BASE="${CATPRED_PIPELINE}/capsule_data_update"
  [[ -d "${CHECKPOINT_BASE}/production" ]] && CHECKPOINT_BASE="${CHECKPOINT_BASE}/production"
fi

# CatPred conda env
if ! $COMMAND_PKG run -n catpred python -c "import sys; sys.exit(0)" 2>/dev/null; then
  echo "Creating conda environment: catpred"
  $COMMAND_PKG env create -f "${CATPRED_REPO}/environment.yml" -n catpred -y
fi
echo "Installing CatPred package in env catpred..."
(cd "${CATPRED_REPO}" && $COMMAND_PKG run -n catpred pip install -e . -q)

# Env file for BEES
ENV_FILE="${REPO_ROOT}/.env.bees"
cat > "${ENV_FILE}" << EOF
# Source before running BEES with kinetics_estimator: catpred
#   source .env.bees
export CATPRED_DIR="${CATPRED_REPO}"
export CATPRED_CHECKPOINT_BASE="${CHECKPOINT_BASE}"
export CATPRED_CONDA_ENV=catpred
EOF
echo "Wrote ${ENV_FILE}"
echo "CatPred install complete. CATPRED_CHECKPOINT_BASE=${CHECKPOINT_BASE}"
