# BEES Makefile (mirrors ARC: https://github.com/ReactionMechanismGenerator/ARC)
# See also: https://reactionmechanismgenerator.github.io/ARC/installation.html

DEVTOOLS_DIR := devtools

.PHONY: all help install install-all install-bees install-catpred test check-env

all: help

help:
	@echo "BEES - Biochemical Engine for Enzymatic kinetic modelS"
	@echo ""
	@echo "Installation (run from BEES repo root):"
	@echo "  make install        Install BEES + CatPred (one-click, same as ./install.sh)"
	@echo "  make install-all   Same as install"
	@echo "  make install-bees  Create/update bees_env only (no CatPred)"
	@echo "  make install-catpred  Install CatPred only (clone, data, catpred env, .env.bees)"
	@echo ""
	@echo "After install:"
	@echo "  conda activate bees_env"
	@echo "  (BEES auto-loads .env.bees when present)"
	@echo "  python BEES.py -i projects/minimal/input.yml"
	@echo ""
	@echo "Other:"
	@echo "  make test          Run pytest"
	@echo "  make check-env     Show Python/conda info"

install: install-all

install-all:
	@echo "Installing BEES and CatPred..."
	bash $(DEVTOOLS_DIR)/install_all.sh

install-bees:
	@echo "Installing BEES only..."
	bash $(DEVTOOLS_DIR)/install_all.sh --no-catpred

install-catpred:
	@echo "Installing CatPred..."
	bash $(DEVTOOLS_DIR)/install_catpred.sh

test:
	conda run -n bees_env python -m pytest tests/ -v

check-env:
	@echo "Python: $$(which python 2>/dev/null || true)"
	@echo "Conda envs:"; conda env list 2>/dev/null || true
	@echo "CATPRED_DIR: $${CATPRED_DIR:-not set}"
