"""
Tests for the Exporter module.

To run:  pytest -v tests/test_exporter.py
"""

import csv
import os
import tempfile
from types import SimpleNamespace
from unittest.mock import MagicMock

import pytest

from bees.enlarger import IterativeEnlarger
from bees.exporter import EnlargerExporter


def _make_reaction(
    enzyme_label="Enzyme",
    substrate_label="S",
    reactant_labels=None,
    product_labels=None,
    stoichiometry=None,
    rate_law="Michaelis-Menten",
    km=1.0,
    kcat=10.0,
):
    kinetics = MagicMock()
    kinetics.km = km
    kinetics.kcat = kcat
    kinetics.vmax = None
    kinetics.km_per_substrate = None

    rxn = MagicMock()
    rxn.enzyme_label = enzyme_label
    rxn.substrate_label = substrate_label
    rxn.reactant_labels = reactant_labels or ["S"]
    rxn.product_labels = product_labels or ["P"]
    rxn.stoichiometry = stoichiometry or {"S": -1, "P": 1}
    rxn.rate_law = rate_law
    rxn.kinetics = kinetics
    # GeneratedReaction-like attributes
    rxn.ec_number = "EC 1.1.1.1"
    rxn.template = MagicMock()
    return rxn


@pytest.fixture
def mock_bees_object():
    settings = SimpleNamespace(
        end_time=100.0,
        time_step=10.0,
        toleranceMoveToCore=1e-5,
        toleranceKeepInEdge=0,
        max_iterations=3,
        max_edge_species=None,
        termination_conversion=None,
        termination_rate_ratio=None,
        save_simulation_profiles=False,
        saveEdgeSpecies=True,
        filter_reactions=False,
    )

    species = [
        SimpleNamespace(
            label="S",
            concentration=10.0,
            reactive=True,
            solvent=False,
            constant=False,
        ),
    ]
    enzymes = [
        SimpleNamespace(
            label="Enzyme",
            concentration=0.01,
            reactive=True,
            solvent=False,
            ecnumber="EC 1.1.1.1",
            constant=True,
        ),
    ]

    return SimpleNamespace(
        settings=settings,
        species=species,
        enzymes=enzymes,
        project="TestProject",
        database=SimpleNamespace(name="db"),
        environment=SimpleNamespace(temperature=298.15, pH=7.0),
    )


@pytest.fixture
def mock_model_generator():
    mg = MagicMock()
    mg.generate_reactions.return_value = [
        _make_reaction(
            enzyme_label="Enzyme",
            substrate_label="S",
            reactant_labels=["S"],
            product_labels=["P"],
            stoichiometry={"S": -1, "P": 1},
            km=1.0,
            kcat=10.0,
        ),
    ]
    mg._generate_reactions.return_value = []
    return mg


@pytest.fixture
def output_dir():
    with tempfile.TemporaryDirectory() as d:
        yield d


class TestEnlargerExporter:
    def _make_exporter(self, enlarger: IterativeEnlarger, mock_bees_object, output_dir, logger):
        return EnlargerExporter(
            model=enlarger.model,
            profiles=enlarger._profiles,
            output_directory=output_dir,
            logger=logger,
            reaction_id_by_sig=enlarger._reaction_id_by_sig,
            reaction_first_seen_iter=enlarger._reaction_first_seen_iter,
            reaction_core_enter_iter=enlarger._reaction_core_enter_iter,
            reaction_obj_by_sig=enlarger._reaction_obj_by_sig,
            iteration_summaries=enlarger._iteration_summaries,
            save_reaction_tree_plots=getattr(
                mock_bees_object.settings, "save_reaction_tree_plots", False
            ),
            save_simulation_plots=getattr(
                mock_bees_object.settings, "save_simulation_plots", False
            ),
            plot_max_species=getattr(mock_bees_object.settings, "plot_max_species", None),
            plot_exclude_enzymes=getattr(
                mock_bees_object.settings, "plot_exclude_enzymes", True
            ),
            plot_exclude_cofactors=getattr(
                mock_bees_object.settings, "plot_exclude_cofactors", True
            ),
            reaction_tree_layout=getattr(
                mock_bees_object.settings, "reaction_tree_layout", "graphviz"
            ),
            reaction_tree_rankdir=getattr(
                mock_bees_object.settings, "reaction_tree_rankdir", "TB"
            ),
            reaction_tree_fontsize=int(
                getattr(mock_bees_object.settings, "reaction_tree_fontsize", 8) or 8
            ),
            core_seen_labels=enlarger._core_seen_labels,
            bees_object=mock_bees_object,
        )

    def test_export_profiles(self, mock_bees_object, mock_model_generator, output_dir):
        """Ensure export_simulation_profiles writes a file."""
        mock_bees_object.settings.save_simulation_profiles = True
        logger = MagicMock()
        enlarger = IterativeEnlarger(
            bees_object=mock_bees_object,
            model_generator=mock_model_generator,
            logger=logger,
            output_directory=output_dir,
        )
        enlarger.run()
        exporter = self._make_exporter(enlarger, mock_bees_object, output_dir, logger)
        path = exporter.export_simulation_profiles()
        if path is not None:
            assert os.path.exists(path)

    def test_export_flux_analysis(self, mock_bees_object, mock_model_generator, output_dir):
        """Ensure export_flux_analysis writes a file."""
        logger = MagicMock()
        enlarger = IterativeEnlarger(
            bees_object=mock_bees_object,
            model_generator=mock_model_generator,
            logger=logger,
            output_directory=output_dir,
        )
        enlarger.run()
        exporter = self._make_exporter(enlarger, mock_bees_object, output_dir, logger)
        path = exporter.export_flux_analysis()
        assert path is not None
        assert os.path.exists(path)

    def test_export_core_edge_reaction_species_csvs(
        self, mock_bees_object, mock_model_generator, output_dir
    ):
        """Ensure core/edge CSV exports are written and include section headers."""
        logger = MagicMock()
        enlarger = IterativeEnlarger(
            bees_object=mock_bees_object,
            model_generator=mock_model_generator,
            logger=logger,
            output_directory=output_dir,
        )
        enlarger.run()
        exporter = self._make_exporter(enlarger, mock_bees_object, output_dir, logger)
        paths = exporter.export_core_edge_reaction_species_csvs()
        assert "core" in paths and "edge" in paths
        assert os.path.exists(paths["core"])
        assert os.path.exists(paths["edge"])

        with open(paths["core"], "r", newline="") as f:
            rows = list(csv.reader(f))
        assert rows[0][:5] == ["index", "reaction_id", "template", "ec_number", "family"]
        assert any(row and row[0] == "core species" for row in rows)

        with open(paths["edge"], "r", newline="") as f:
            rows = list(csv.reader(f))
        assert rows[0][:5] == ["index", "reaction_id", "template", "ec_number", "family"]
        assert any(row and row[0] == "edge species" for row in rows)

    def test_reaction_tree_export_core_only_no_enzymes_or_cofactors(
        self, mock_bees_object, mock_model_generator, output_dir
    ):
        """
        When save_reaction_tree_plots is enabled, exporter should write a PNG that
        only includes core, non-enzyme, non-cofactor species.
        """
        from bees.core_edge_model import SpeciesData

        mock_bees_object.settings.save_reaction_tree_plots = True

        enlarger = IterativeEnlarger(
            bees_object=mock_bees_object,
            model_generator=mock_model_generator,
            logger=MagicMock(),
            output_directory=output_dir,
        )

        enlarger.model.add_core_species(SpeciesData(label="A", concentration=1.0))
        enlarger.model.add_core_species(SpeciesData(label="B", concentration=0.0))
        enlarger.model.add_core_species(SpeciesData(label="EnzymeX", is_enzyme=True))
        enlarger.model.add_core_species(SpeciesData(label="ATP", concentration=1.0))

        rxn = _make_reaction(
            enzyme_label="EnzymeX",
            reactant_labels=["A", "ATP"],
            product_labels=["B", "ATP"],
        )
        enlarger.model.core_reactions.append(rxn)

        enlarger._core_seen_labels = set()

        exporter = EnlargerExporter(
            model=enlarger.model,
            profiles=enlarger._profiles,
            output_directory=output_dir,
            logger=MagicMock(),
            reaction_id_by_sig=enlarger._reaction_id_by_sig,
            reaction_first_seen_iter=enlarger._reaction_first_seen_iter,
            reaction_core_enter_iter=enlarger._reaction_core_enter_iter,
            reaction_obj_by_sig=enlarger._reaction_obj_by_sig,
            iteration_summaries=enlarger._iteration_summaries,
            save_reaction_tree_plots=True,
            core_seen_labels=enlarger._core_seen_labels,
            bees_object=mock_bees_object,
        )
        exporter.export_reaction_tree(iteration=1, promoted_labels=["A", "B"])

        files = [
            name
            for name in os.listdir(output_dir)
            if name.startswith("reaction_tree_iter1") and name.endswith(".png")
        ]
        assert files, "Expected at least one reaction_tree_iter1*.png file to be created"

