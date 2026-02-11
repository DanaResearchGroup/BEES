#!/usr/bin/env bash
# One-click install for BEES and CatPred (same as: make install)
# Run from BEES repo root: ./install.sh


set -e
SCRIPT_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" &>/dev/null && pwd)"
bash "$SCRIPT_DIR/devtools/install_all.sh"

echo ""
echo "To run BEES: conda activate bees_env"
echo "With CatPred: source .env.bees   then run BEES as above."
