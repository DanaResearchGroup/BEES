"""
Tests for the Exporter module.

To run:  pytest -v tests/test_exporter.py
"""

import csv
import os
import tempfile
from types import SimpleNamespace
from unittest.mock import MagicMock

import numpy as np
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
def mock_reaction_generator():
    mg = MagicMock()
    mg._generate_reactions.return_value = [
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
            bees_object=mock_bees_object,
        )

    def test_export_profiles(self, mock_bees_object, mock_reaction_generator, output_dir):
        """Ensure export_simulation_profiles writes a file."""
        mock_bees_object.settings.save_simulation_profiles = True
        logger = MagicMock()
        enlarger = IterativeEnlarger(
            bees_object=mock_bees_object,
            reaction_generator=mock_reaction_generator,
            logger=logger,
            output_directory=output_dir,
        )
        enlarger.run()
        exporter = self._make_exporter(enlarger, mock_bees_object, output_dir, logger)
        path = exporter.export_simulation_profiles()
        assert path is not None and os.path.exists(path)

    def test_export_flux_analysis(self, mock_bees_object, mock_reaction_generator, output_dir):
        """Ensure export_flux_analysis writes a file."""
        logger = MagicMock()
        enlarger = IterativeEnlarger(
            bees_object=mock_bees_object,
            reaction_generator=mock_reaction_generator,
            logger=logger,
            output_directory=output_dir,
        )
        enlarger.run()
        exporter = self._make_exporter(enlarger, mock_bees_object, output_dir, logger)
        path = exporter.export_flux_analysis()
        assert path is not None
        assert os.path.exists(path)

    def test_export_core_edge_reaction_species_csvs(
        self, mock_bees_object, mock_reaction_generator, output_dir
    ):
        """Ensure core/edge CSV exports are written and include section headers."""
        logger = MagicMock()
        enlarger = IterativeEnlarger(
            bees_object=mock_bees_object,
            reaction_generator=mock_reaction_generator,
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
        assert rows[1][1] == "R1" and rows[1][-2:] == ["S", "P"]
        assert any(row and row[0] == "core species" for row in rows)

        with open(paths["edge"], "r", newline="") as f:
            rows = list(csv.reader(f))
        assert rows[0][:5] == ["index", "reaction_id", "template", "ec_number", "family"]
        assert any(row and row[0] == "edge species" for row in rows)

    def test_edge_species_csv_has_peak_rate_ratio(self, output_dir):
        """Edge species section carries the peak |R|/R_char and its iteration (blank if unknown)."""
        edge = [
            SimpleNamespace(label="Butanoyl-CoA", is_enzyme=False, constant=False),
            SimpleNamespace(label="acetate", is_enzyme=False, constant=False),
        ]
        model = SimpleNamespace(core_reactions=[], edge_reactions=[], core_species=[], edge_species=edge)
        exporter = EnlargerExporter(
            model=model,
            profiles=[],
            output_directory=output_dir,
            logger=MagicMock(),
            reaction_id_by_sig={},
            reaction_first_seen_iter={},
            reaction_core_enter_iter={},
            reaction_obj_by_sig={},
            iteration_summaries=[],
            edge_peak_rate_ratio={"butanoyl-coa": (3.1e-3, 42)},
        )
        paths = exporter.export_core_edge_reaction_species_csvs()

        with open(paths["edge"], "r", newline="") as f:
            rows = list(csv.reader(f))
        start = next(i for i, row in enumerate(rows) if row and row[0] == "edge species")
        assert rows[start + 1] == ["index", "label", "is_enzyme", "constant", "peak_rate_ratio", "peak_iteration"]
        assert rows[start + 2] == ["1", "Butanoyl-CoA", "False", "False", "3.100000e-03", "42"]
        assert rows[start + 3] == ["2", "acetate", "False", "False", "", ""]

        with open(paths["core"], "r", newline="") as f:
            rows = list(csv.reader(f))
        assert ["index", "label", "is_enzyme", "constant"] in rows

# ---------------------------------------------------------------------------
# COPASI-compatibility: SBML `reversible` must match the kinetic law emitted
# (fallback-thermo reactions must be reversible="false").
# ---------------------------------------------------------------------------

def _make_reaction_with_thermo(
    *,
    enzyme_label="EnzymeR",
    reactant_labels=("A",),
    product_labels=("B",),
    stoichiometry=None,
    km=1.0,
    kcat=10.0,
    thermo_source="equilibrator",
    irreversible=False,
    keq=100.0,
    kcat_rev=0.1,
    dgr=-10.0,
    template_reversible=True,
):
    """Build a reaction with a real-enough .thermo attribute for SBML export."""
    rxn = _make_reaction(
        enzyme_label=enzyme_label,
        reactant_labels=list(reactant_labels),
        product_labels=list(product_labels),
        stoichiometry=stoichiometry
        or {**{r: -1 for r in reactant_labels}, **{p: 1 for p in product_labels}},
        km=km,
        kcat=kcat,
    )
    thermo = SimpleNamespace(
        source=thermo_source,
        irreversible=irreversible,
        keq=keq,
        kcat_rev=kcat_rev,
        dgr_prime_kJmol=dgr,
        sigma_kJmol=1.0,
    )
    rxn.thermo = thermo
    rxn.template.reversible = template_reversible
    return rxn


def _export_sbml_with_reactions(
    rxns, output_dir, mock_bees_object, *, extra_species=(), strict_invariant=True,
    exporter_logger=None, add_reaction_species=True, explored=(),
):
    """Build an enlarger, attach reactions and species, export SBML, return path.

    ``exporter_logger`` overrides the MagicMock logger on the exporter — pass a
    strict stub to assert the export uses single-arg log calls (BEES's real
    ``Logger`` signature), which a MagicMock would silently swallow.
    """
    from bees.core_edge_model import SpeciesData

    mg = MagicMock()
    mg.generate_reactions.return_value = []
    mg._generate_reactions.return_value = []

    enlarger = IterativeEnlarger(
        bees_object=mock_bees_object,
        reaction_generator=mg,
        logger=MagicMock(),
        output_directory=output_dir,
    )

    # Add every species referenced by any reaction (plus the enzyme).
    species_seen = set()
    for rxn in rxns:
        for lab in list(rxn.reactant_labels) + list(rxn.product_labels):
            if not add_reaction_species or lab in species_seen:
                continue
            species_seen.add(lab)
            enlarger.model.add_core_species(SpeciesData(label=lab, concentration=1.0))
        enz_lab = rxn.enzyme_label
        if enz_lab not in species_seen:
            species_seen.add(enz_lab)
            enlarger.model.add_core_species(
                SpeciesData(label=enz_lab, concentration=0.01, is_enzyme=True)
            )
        enlarger.model.core_reactions.append(rxn)
    for sd in extra_species:
        enlarger.model.add_core_species(sd)
    enlarger.model.explored_species_lc.update(lab.lower() for lab in explored)


    exporter = EnlargerExporter(
        model=enlarger.model,
        profiles=enlarger._profiles,
        output_directory=output_dir,
        logger=exporter_logger if exporter_logger is not None else MagicMock(),
        reaction_id_by_sig=enlarger._reaction_id_by_sig,
        reaction_first_seen_iter=enlarger._reaction_first_seen_iter,
        reaction_core_enter_iter=enlarger._reaction_core_enter_iter,
        reaction_obj_by_sig=enlarger._reaction_obj_by_sig,
        iteration_summaries=enlarger._iteration_summaries,
        bees_object=mock_bees_object,
    )
    path = exporter.export_sbml(filename="model.xml", strict_invariant=strict_invariant)
    return path


def test_fallback_thermo_reaction_is_irreversible_in_sbml(
    mock_bees_object, output_dir
):
    """
    A reaction whose thermo is `source="fallback"`
    (eQuilibrator could not resolve it) must be exported with
    `reversible="false"`, because the kinetic law emitted is forward-only MM.
    SBML attribute and kinetic-law form must agree, or COPASI raises
    "never negative — unexpected for reversible reaction".
    """
    libsbml = pytest.importorskip("libsbml")

    rxn_fallback = _make_reaction_with_thermo(
        enzyme_label="EnzymeF",
        reactant_labels=("X",),
        product_labels=("Y",),
        thermo_source="fallback",
        irreversible=True,
        keq=float("nan"),
        kcat_rev=None,
        dgr=float("nan"),
    )

    # Toy 1-reaction fixture: Y has no consumer. The export invariant is
    # designed to catch this in real models; relax it here so the test can
    # focus on its actual concern (reversibility flag on fallback thermo).
    path = _export_sbml_with_reactions(
        [rxn_fallback], output_dir, mock_bees_object, strict_invariant=False,
    )
    assert path is not None and os.path.exists(path)

    doc = libsbml.SBMLReader().readSBMLFromFile(path)
    m = doc.getModel()
    assert m.getNumReactions() == 1
    sbml_rxn = m.getReaction(0)
    assert sbml_rxn.getReversible() is False, (
        "Fallback-thermo reaction must be reversible=false in SBML to match "
        "its forward-only kinetic law"
    )


def test_kinetic_law_is_scaled_by_compartment_volume(mock_bees_object, output_dir):
    """
    SBML defines a kineticLaw as a *substance* rate (amount/time) and a solver
    computes dC/dt = kineticLaw / V. BEES builds *concentration*-rate formulas
    (mM/s), so each kineticLaw must be multiplied by the compartment volume to
    survive that division — mirroring RMG's liquid solver (`core_species_rates *
    V`). Without it the exported model integrates 1/V times too fast in COPASI.

    Guards against regressing to a bare concentration-rate formula or a
    compartment size that does not cancel.
    """
    libsbml = pytest.importorskip("libsbml")

    rxn = _make_reaction_with_thermo(
        enzyme_label="EnzymeV",
        reactant_labels=("A",),
        product_labels=("B",),
    )
    path = _export_sbml_with_reactions([rxn], output_dir, mock_bees_object)
    assert path is not None and os.path.exists(path)

    doc = libsbml.SBMLReader().readSBMLFromFile(path)
    m = doc.getModel()

    comp = m.getCompartment("compartment1")
    assert comp is not None and comp.getSize() == 1.0, (
        "compartment 'compartment1' must have size 1.0 so amount (mmol) == "
        "concentration (mM) and dC/dt equals the bare rate expression"
    )

    assert m.getNumReactions() >= 1
    for i in range(m.getNumReactions()):
        kl = m.getReaction(i).getKineticLaw()
        assert kl is not None, "every reaction should have a kineticLaw"
        ci_names = {n for n in _walk_ast_ci_names(kl.getMath())}
        assert "compartment1" in ci_names, (
            "kineticLaw must reference the compartment 'compartment1' as a volume "
            "factor; otherwise the exported model runs 1/V too fast in COPASI"
        )


def _walk_ast_ci_names(ast_node):
    """Yield the names of every <ci> (identifier) node in an AST subtree."""
    if ast_node is None:
        return
    if ast_node.isName() and ast_node.getName():
        yield ast_node.getName()
    for i in range(ast_node.getNumChildren()):
        yield from _walk_ast_ci_names(ast_node.getChild(i))


def test_feedback_inhibition_exported_to_kinetic_law(mock_bees_object, output_dir):
    """
    A reaction carrying `feedback_inhibitors` (Enzyme.feedback_inhibition, the
    simulator's shared-site factor) must export the multiplier 1/(1+Σ([I]/Ki)^h) into its
    kineticLaw, declare each inhibitor as a modifier, and emit Ki/hill parameters.
    Without this the SBML shows the baseline shape, not the feedback shape.
    """
    libsbml = pytest.importorskip("libsbml")
    from bees.core_edge_model import SpeciesData

    rxn = _make_reaction_with_thermo(
        enzyme_label="FabH",
        reactant_labels=("A",),
        product_labels=("B",),
        thermo_source="fallback",
        irreversible=True,
        keq=float("nan"),
        kcat_rev=None,
        dgr=float("nan"),
    )
    rxn.feedback_inhibitors = {"inhibitor_x": (0.005, 1.0)}

    path = _export_sbml_with_reactions(
        [rxn], output_dir, mock_bees_object, strict_invariant=False,
        extra_species=(SpeciesData(label="inhibitor_x", concentration=0.0),),
    )
    doc = libsbml.SBMLReader().readSBMLFromFile(path)
    m = doc.getModel()
    rx = m.getReaction(0)

    mods = {rx.getModifier(i).getSpecies() for i in range(rx.getNumModifiers())}
    assert "inhibitor_x" in mods, "inhibitor must be a modifierSpeciesReference"

    ki = next(
        (m.getParameter(i) for i in range(m.getNumParameters())
         if m.getParameter(i).getId().startswith("Ki_fb_")),
        None,
    )
    assert ki is not None and ki.getValue() == 0.005, "Ki parameter must be emitted"

    ci_names = set(_walk_ast_ci_names(rx.getKineticLaw().getMath()))
    assert "inhibitor_x" in ci_names, "kineticLaw must reference the inhibitor"
    assert any(n.startswith("Ki_fb_") for n in ci_names), (
        "kineticLaw must reference the feedback Ki parameter"
    )


def test_sbml_passes_libsbml_consistency_check(mock_bees_object, output_dir):
    """The exported document has no SBML errors (e.g. unit ids clashing with built-in units)."""
    libsbml = pytest.importorskip("libsbml")
    from bees.core_edge_model import SpeciesData

    rxn = _make_reaction_with_thermo(
        enzyme_label="FabH", reactant_labels=("A",), product_labels=("B",),
        thermo_source="fallback", irreversible=True, keq=float("nan"),
        kcat_rev=None, dgr=float("nan"),
    )
    rxn.feedback_inhibitors = {"inh_x": (0.005, 1.0)}
    path = _export_sbml_with_reactions(
        [rxn], output_dir, mock_bees_object, strict_invariant=False,
        extra_species=(SpeciesData(label="inh_x", concentration=0.0),),
    )

    doc = libsbml.SBMLReader().readSBMLFromFile(path)
    doc.setConsistencyChecks(libsbml.LIBSBML_CAT_UNITS_CONSISTENCY, False)
    doc.checkConsistency()
    errors = [
        doc.getError(i).getMessage().strip()
        for i in range(doc.getNumErrors())
        if doc.getError(i).getSeverity() >= libsbml.LIBSBML_SEV_ERROR
    ]
    assert errors == []
    m = doc.getModel()
    assert (m.getTimeUnits(), m.getVolumeUnits(), m.getSubstanceUnits()) == ("second", "litre", "mmol")


def test_feedback_inhibitors_exported_as_shared_site(mock_bees_object, output_dir):
    """Two feedback inhibitors export as one shared site, matching the simulator."""
    libsbml = pytest.importorskip("libsbml")
    from bees.core_edge_model import SpeciesData
    from bees.flux_calculator import compute_mm_rate

    rxn = _make_reaction_with_thermo(
        enzyme_label="FabH", reactant_labels=("A",), product_labels=("B",),
        thermo_source="fallback", irreversible=True, keq=float("nan"),
        kcat_rev=None, dgr=float("nan"),
    )
    rxn.feedback_inhibitors = {"inh_x": (0.005, 1.0), "inh_y": (0.005, 1.0)}

    path = _export_sbml_with_reactions(
        [rxn], output_dir, mock_bees_object, strict_invariant=False,
        extra_species=(
            SpeciesData(label="inh_x", concentration=0.0),
            SpeciesData(label="inh_y", concentration=0.0),
        ),
    )
    doc = libsbml.SBMLReader().readSBMLFromFile(path)
    m = doc.getModel()
    conc = {"a": 0.1, "b": 0.0, "inh_x": 0.005, "inh_y": 0.005}
    ns = {p.getId(): p.getValue() for p in m.getListOfParameters()}
    ns["compartment1"] = 1.0
    for s in m.getListOfSpecies():
        ns[s.getId()] = conc.get(s.getName().lower(), 0.01)
    base = compute_mm_rate(rxn, conc, {"fabh": ns.get(next(
        s.getId() for s in m.getListOfSpecies() if s.getName().lower() == "fabh"
    ), 0.01)})

    formula = libsbml.formulaToL3String(m.getReaction(0).getKineticLaw().getMath())
    v = eval(formula.replace("^", "**"), {"__builtins__": {}}, ns)
    assert v == pytest.approx(base / 3.0, rel=1e-9), formula


def test_sbml_kinetic_laws_evaluate_like_simulator(mock_bees_object, output_dir):
    """
    Each exported kineticLaw must give the simulator's rate. The reversible reaction
    carries a scaled substrate Km and a kcat_rev that is not Haldane-consistent with it
    (as after hydrophobic_chain_length_km); the simulator uses 1 - Q/Keq, so must the SBML.
    """
    libsbml = pytest.importorskip("libsbml")
    from bees.flux_calculator import compute_mm_rate, compute_reversible_mm_rate

    rxn_rev = _make_reaction_with_thermo(
        enzyme_label="EnzymeR", reactant_labels=("A",), product_labels=("B", "H2O"),
        stoichiometry={"A": -1, "B": 1, "H2O": 1}, keq=1.19, kcat=6.5, kcat_rev=4.67,
    )
    rxn_rev.kinetics.km_per_substrate = {"A": 1.4e-5, "B": 0.039}
    rxn_pi = _make_reaction_with_thermo(
        enzyme_label="EnzymeP", reactant_labels=("B",), product_labels=("C",),
        stoichiometry={"B": -1, "C": 1}, irreversible=True, kcat=2.0,
    )
    rxn_pi.kinetics.km_per_substrate = {"B": 0.01, "C": 0.05}

    path = _export_sbml_with_reactions(
        [rxn_rev, rxn_pi], output_dir, mock_bees_object, strict_invariant=False,
    )
    doc = libsbml.SBMLReader().readSBMLFromFile(path)
    m = doc.getModel()
    conc = {"a": 1.9e-6, "b": 2.27e-6, "h2o": 55.0, "c": 0.02}
    enz = {"enzymer": 0.01, "enzymep": 0.01}
    ns = {p.getId(): p.getValue() for p in m.getListOfParameters()}
    ns["compartment1"] = 1.0
    for s in m.getListOfSpecies():
        ns[s.getId()] = conc.get(s.getName().lower(), 0.01)

    expected = {
        "EnzymeR": compute_reversible_mm_rate(rxn_rev, conc, enz),
        "EnzymeP": compute_mm_rate(rxn_pi, conc, enz),
    }
    assert expected["EnzymeR"] < 1e-3 * 6.5 * 0.01  # near equilibrium, unlike kcat_rev form
    for i in range(m.getNumReactions()):
        rx = m.getReaction(i)
        enzyme = rx.getName().split(":")[0]
        formula = libsbml.formulaToL3String(rx.getKineticLaw().getMath()).replace("^", "**")
        v = eval(formula, {"__builtins__": {}}, ns)
        assert v == pytest.approx(expected[enzyme], rel=1e-9), (enzyme, formula)
    assert m.getReaction(0).getReversible() is True
    assert m.getReaction(1).getReversible() is False


def test_sbml_reversible_law_uses_1_molar_standard_state(mock_bees_object, output_dir):
    """A + B <=> C with Keq = 10 (1 M standard state) is at equilibrium at
    A = B = 1 mM, C = 0.01 mM. The SBML law must vanish there and match the
    simulator away from it."""
    libsbml = pytest.importorskip("libsbml")
    from bees.flux_calculator import compute_reversible_mm_rate

    rxn = _make_reaction_with_thermo(
        enzyme_label="EnzymeR", reactant_labels=("A", "B"), product_labels=("C",),
        stoichiometry={"A": -1, "B": -1, "C": 1}, keq=10.0, kcat=6.5, kcat_rev=None,
    )
    rxn.kinetics.km_per_substrate = {"A": 0.1, "B": 0.1, "C": 0.1}
    path = _export_sbml_with_reactions([rxn], output_dir, mock_bees_object, strict_invariant=False)
    doc = libsbml.SBMLReader().readSBMLFromFile(path)
    m = doc.getModel()
    formula = libsbml.formulaToL3String(m.getReaction(0).getKineticLaw().getMath()).replace("^", "**")
    enz = {"enzymer": 0.01}
    for conc in ({"a": 1.0, "b": 1.0, "c": 0.01}, {"a": 1.0, "b": 1.0, "c": 0.02},
                 {"a": 0.3, "b": 2.0, "c": 0.0}):
        ns = {p.getId(): p.getValue() for p in m.getListOfParameters()}
        ns["compartment1"] = 1.0
        for s in m.getListOfSpecies():
            ns[s.getId()] = conc.get(s.getName().lower(), 0.01)
        v = eval(formula, {"__builtins__": {}}, ns)
        expected = compute_reversible_mm_rate(rxn, conc, enz)
        assert v == pytest.approx(expected, rel=1e-9, abs=1e-15), (conc, formula)
    at_eq = {"a": 1.0, "b": 1.0, "c": 0.01}
    assert compute_reversible_mm_rate(rxn, at_eq, enz) == 0.0


class _StrictLogger:
    """Mimics bees.logger.Logger: every level takes a single message string.

    A MagicMock accepts `logger.debug(msg, a, b)` silently; the real BEES logger
    raises TypeError. This stub reproduces that so tests catch %-style log calls.
    """

    def _check(self, message, *args):
        assert not args, (
            "BEES logger takes a single pre-formatted message; use an f-string, "
            f"not %-args. Got extra args: {args!r}"
        )

    debug = info = warning = error = _check


def test_export_uses_single_arg_log_calls_on_missing_inhibitor(
    mock_bees_object, output_dir
):
    """
    Regression: the feedback-export path logs when an inhibitor is not a model
    species. That log call (and all others on the export path) must use BEES's
    single-arg Logger signature, not %-style args. Exercised with a strict logger
    and a feedback inhibitor that is deliberately absent from the species set.
    """
    pytest.importorskip("libsbml")

    rxn = _make_reaction_with_thermo(
        enzyme_label="FabH",
        reactant_labels=("A",),
        product_labels=("B",),
        thermo_source="fallback",
        irreversible=True,
        keq=float("nan"),
        kcat_rev=None,
        dgr=float("nan"),
    )
    # Inhibitor "ghost_acid" is never added as a species -> triggers the debug log.
    rxn.feedback_inhibitors = {"ghost_acid": (0.005, 1.0)}

    # Must not raise (a %-style log call would raise via _StrictLogger).
    path = _export_sbml_with_reactions(
        [rxn], output_dir, mock_bees_object, strict_invariant=False,
        exporter_logger=_StrictLogger(),
    )
    assert path is not None and os.path.exists(path)


def test_export_uses_single_arg_log_calls_on_missing_kinetics(
    mock_bees_object, output_dir
):
    """A reaction without kinetics is exported without a rate law, not a crash."""
    libsbml = pytest.importorskip("libsbml")

    rxn = _make_reaction_with_thermo(enzyme_label="E1")
    rxn.kinetics = None
    path = _export_sbml_with_reactions(
        [rxn], output_dir, mock_bees_object, strict_invariant=False,
        exporter_logger=_StrictLogger(),
    )
    assert path is not None and os.path.exists(path)
    doc = libsbml.SBMLReader().readSBMLFromFile(path)
    m = doc.getModel()
    assert m.getNumReactions() == 1
    assert not m.getReaction(0).isSetKineticLaw()


def test_export_resolves_ontology_alias_like_simulator(mock_bees_object, output_dir):
    """A reaction written with an alias (D-glucose) uses the model species (Glucose)."""
    libsbml = pytest.importorskip("libsbml")
    from bees.core_edge_model import SpeciesData

    rxn = _make_reaction_with_thermo(enzyme_label="E1", reactant_labels=("D-glucose",), product_labels=("B",))
    path = _export_sbml_with_reactions(
        [rxn], output_dir, mock_bees_object, add_reaction_species=False,
        extra_species=[SpeciesData(label="Glucose", concentration=1.0), SpeciesData(label="B", concentration=0.0)],
    )
    doc = libsbml.SBMLReader().readSBMLFromFile(path)
    r = doc.getModel().getReaction(0)
    assert [r.getReactant(i).getSpecies() for i in range(r.getNumReactants())] == ["Glucose"]


def test_export_accepts_concentration_range(mock_bees_object, output_dir):
    """A (min, max) concentration exports the lower bound, as the simulator uses."""
    libsbml = pytest.importorskip("libsbml")

    mock_bees_object.species[0].label = "A"
    mock_bees_object.species[0].concentration = (2.0, 5.0)
    mock_bees_object.enzymes[0].label = "EnzymeR"
    mock_bees_object.enzymes[0].concentration = (0.01, 0.02)
    rxn = _make_reaction_with_thermo()
    path = _export_sbml_with_reactions(
        [rxn], output_dir, mock_bees_object, strict_invariant=False,
    )
    doc = libsbml.readSBMLFromFile(path)
    model = doc.getModel()
    a = next(model.getSpecies(i) for i in range(model.getNumSpecies())
             if model.getSpecies(i).getName() == "A")
    assert a.getInitialConcentration() == pytest.approx(2.0)
    enz = next(model.getParameter(i) for i in range(model.getNumParameters())
               if model.getParameter(i).getName() == "[EnzymeR]")
    assert enz.getValue() == pytest.approx(0.01)


def _free_acid_reactions():
    rxn_c8 = _make_reaction_with_thermo(
        enzyme_label="TesA8", reactant_labels=("octanoyl-[ACP]",),
        product_labels=("octanoate",), thermo_source="fallback",
        irreversible=True, keq=float("nan"), kcat_rev=None, dgr=float("nan"),
    )
    rxn_c16 = _make_reaction_with_thermo(
        enzyme_label="TesA16", reactant_labels=("hexadecanoyl-[ACP]",),
        product_labels=("hexadecanoate",), thermo_source="fallback",
        irreversible=True, keq=float("nan"), kcat_rev=None, dgr=float("nan"),
    )
    return [rxn_c8, rxn_c16]


def test_palmitic_equivalents_needs_a_fas_calibration(mock_bees_object, output_dir):
    """A generic run (no FAS calibration enabled) writes no FAS observable."""
    libsbml = pytest.importorskip("libsbml")
    path = _export_sbml_with_reactions(
        _free_acid_reactions(), output_dir, mock_bees_object, strict_invariant=False,
    )
    doc = libsbml.SBMLReader().readSBMLFromFile(path)
    assert doc.getModel().getParameter("palmitic_equivalents_uM") is None


def test_free_fatty_acid_global_quantities_exported(mock_bees_object, output_dir):
    """
    With a FAS calibration enabled, the exporter must emit one global quantity
    (parameter + assignment rule) for the Yu-2011 / Ruppe S2A time-course PE metric:
      palmitic_equivalents_uM = 1000·Σ((carbons-2)/14)·[acid]
    (malonyl-label weight: Yu et al. 2011 counts [2-14C]malonyl, palmitate = 7 units).
    """
    libsbml = pytest.importorskip("libsbml")

    mock_bees_object.settings.calibrations = ["tesa_long_chain_preference"]
    path = _export_sbml_with_reactions(
        _free_acid_reactions(), output_dir, mock_bees_object, strict_invariant=False,
    )
    doc = libsbml.SBMLReader().readSBMLFromFile(path)
    m = doc.getModel()

    p = m.getParameter("palmitic_equivalents_uM")
    assert p is not None and p.getConstant() is False, (
        "palmitic_equivalents_uM must be a non-constant global quantity"
    )
    rule_vars = {m.getRule(i).getVariable() for i in range(m.getNumRules())}
    assert "palmitic_equivalents_uM" in rule_vars

    palm_rule = next(
        m.getRule(i) for i in range(m.getNumRules())
        if m.getRule(i).getVariable() == "palmitic_equivalents_uM"
    )
    ci_names = set(_walk_ast_ci_names(palm_rule.getMath()))
    assert {"octanoate", "hexadecanoate"} <= ci_names, (
        "palmitic-equivalents rule must weight every free-acid species"
    )

    formula = libsbml.formulaToL3String(palm_rule.getMath())

    def _rule_value(octanoate, hexadecanoate):
        return eval(formula, {"__builtins__": {}},
                    {"octanoate": octanoate, "hexadecanoate": hexadecanoate})

    assert _rule_value(1.0, 0.0) == pytest.approx(1000.0 * (8 - 2) / 14), (
        "C8 term must carry the malonyl-label weight (8-2)/14"
    )
    assert _rule_value(0.0, 1.0) == pytest.approx(1000.0 * 14 / 14), (
        "C16 term must carry weight 14/14 (palmitate = 1)"
    )


# ---------------------------------------------------------------------------
# Production-only-species invariant (RMG-style edge isolation companion).
# A non-boundary, non-terminal species that appears as a product in some
# reaction but never as a reactant will accumulate without bound under
# integration. The export must refuse to write such a model by default and
# only warn when explicitly demoted.
# ---------------------------------------------------------------------------

class TestProductionOnlyInvariant:
    def _make_pair(self, *, product_label):
        """A -> product_label, plus a downstream consumer of A so A itself isn't a leak."""
        # consumer keeps A from being production-only
        seed = _make_reaction_with_thermo(
            reactant_labels=("Seed",),
            product_labels=("A",),
            template_reversible=False,
            irreversible=True,
        )
        leaker = _make_reaction_with_thermo(
            enzyme_label="EnzL",
            reactant_labels=("A",),
            product_labels=(product_label,),
            template_reversible=False,
            irreversible=True,
        )
        return [seed, leaker]

    def test_refuses_production_only_species(self, output_dir, mock_bees_object):
        # "Mystery" is produced (by leaker) but never consumed and is not on
        # any allowlist or suffix family → must raise.
        rxns = self._make_pair(product_label="MysterySink")
        with pytest.raises(ValueError, match="produced but never consumed"):
            _export_sbml_with_reactions(rxns, output_dir, mock_bees_object)

    def test_allows_searched_end_products(self, output_dir, mock_bees_object):
        # Searched for consumers and none found (e.g. a free fatty acid) → end product.
        rxns = self._make_pair(product_label="hexanoate")
        path = _export_sbml_with_reactions(
            rxns, output_dir, mock_bees_object, explored=("hexanoate",),
        )
        assert path is not None and os.path.exists(path)

    def test_name_ending_in_ate_is_not_exempt(self, output_dir, mock_bees_object):
        # Never searched: an "-ate" name alone does not hide the leak.
        rxns = self._make_pair(product_label="pyruvate")
        with pytest.raises(ValueError, match="produced but never consumed"):
            _export_sbml_with_reactions(rxns, output_dir, mock_bees_object)

    def test_allows_general_cofactor_terminals(self, output_dir, mock_bees_object):
        # CoA is a general cofactor — production-only cofactors are not leaks.
        rxns = self._make_pair(product_label="Coenzyme A")
        path = _export_sbml_with_reactions(rxns, output_dir, mock_bees_object)
        assert path is not None and os.path.exists(path)

    def test_strict_invariant_false_demotes_to_warning(
        self, output_dir, mock_bees_object,
    ):
        # Same broken model as test_refuses…, but called with strict=False.
        # Build the exporter manually so we can pass strict_invariant=False.
        from bees.core_edge_model import SpeciesData
        from bees.enlarger import IterativeEnlarger

        rxns = self._make_pair(product_label="MysterySink")
        mg = MagicMock()
        mg.generate_reactions.return_value = []
        mg._generate_reactions.return_value = []
        enlarger = IterativeEnlarger(
            bees_object=mock_bees_object,
            reaction_generator=mg,
            logger=MagicMock(),
            output_directory=output_dir,
        )
        species_seen = set()
        for rxn in rxns:
            for lab in list(rxn.reactant_labels) + list(rxn.product_labels):
                if lab in species_seen:
                    continue
                species_seen.add(lab)
                enlarger.model.add_core_species(SpeciesData(label=lab, concentration=1.0))
            enz_lab = rxn.enzyme_label
            if enz_lab not in species_seen:
                species_seen.add(enz_lab)
                enlarger.model.add_core_species(
                    SpeciesData(label=enz_lab, concentration=0.01, is_enzyme=True)
                )
            enlarger.model.core_reactions.append(rxn)

        logger = MagicMock()
        exporter = EnlargerExporter(
            model=enlarger.model,
            profiles=enlarger._profiles,
            output_directory=output_dir,
            logger=logger,
            reaction_id_by_sig=enlarger._reaction_id_by_sig,
            reaction_first_seen_iter=enlarger._reaction_first_seen_iter,
            reaction_core_enter_iter=enlarger._reaction_core_enter_iter,
            reaction_obj_by_sig=enlarger._reaction_obj_by_sig,
            iteration_summaries=enlarger._iteration_summaries,
            bees_object=mock_bees_object,
        )
        path = exporter.export_sbml(filename="model.xml", strict_invariant=False)
        assert path is not None and os.path.exists(path)
        # The warning should have been logged.
        warning_msgs = [str(c) for c in logger.warning.call_args_list]
        assert any("produced but never consumed" in m for m in warning_msgs), (
            f"expected production-only warning, got: {warning_msgs}"
        )


# ---------------------------------------------------------------------------
# Inhibition and dosing in the SBML export, integrated independently of BEES
# ---------------------------------------------------------------------------

def _integrate_sbml(path, t_eval, rtol=1e-10, atol=1e-13):
    """Integrate an exported SBML file from its own math (libsbml AST -> Python), not from BEES.

    Handles what BEES writes: mass-action stoichiometry, kinetic laws over species and
    parameters, the `time` symbol, piecewise/&&, and `time >= t` events that add to a species.
    Returns {species name: concentrations at t_eval} (compartment size 1).
    """
    import math
    import re

    from scipy.integrate import solve_ivp

    libsbml = pytest.importorskip("libsbml")
    doc = libsbml.readSBML(path)
    assert doc.getNumErrors(libsbml.LIBSBML_SEV_ERROR) == 0, doc.getErrorLog().toString()
    model = doc.getModel()

    ids = [model.getSpecies(i).getId() for i in range(model.getNumSpecies())]
    pos = {sid: i for i, sid in enumerate(ids)}
    fixed = {sid for sid in ids if model.getSpecies(sid).getBoundaryCondition()}
    y0 = np.array([model.getSpecies(sid).getInitialConcentration() for sid in ids])
    consts = {model.getParameter(i).getId(): model.getParameter(i).getValue() for i in range(model.getNumParameters())}
    consts["compartment1"] = 1.0

    def pw(*a):
        for val, cond in zip(a[0::2], a[1::2]):
            if cond:
                return val
        return a[-1]

    def compile_expr(ast):
        text = libsbml.formulaToL3String(ast).replace("^", "**").replace("&&", " and ")
        text = re.sub(r"\bpiecewise\(", "_pw(", text)
        return compile(text, "<sbml>", "eval"), text

    reactions = []
    breaks = set()
    for r in (model.getReaction(i) for i in range(model.getNumReactions())):
        code, text = compile_expr(r.getKineticLaw().getMath())
        breaks.update(float(x) for x in re.findall(r"time (?:>=|<) ([0-9.eE+-]+)", text))
        stoich = {}
        for k in range(r.getNumReactants()):
            ref = r.getReactant(k)
            stoich[ref.getSpecies()] = stoich.get(ref.getSpecies(), 0.0) - ref.getStoichiometry()
        for k in range(r.getNumProducts()):
            ref = r.getProduct(k)
            stoich[ref.getSpecies()] = stoich.get(ref.getSpecies(), 0.0) + ref.getStoichiometry()
        reactions.append((code, stoich))

    events = []
    for ev in (model.getEvent(i) for i in range(model.getNumEvents())):
        t_ev = float(re.search(r"time >= ([0-9.eE+-]+)", libsbml.formulaToL3String(ev.getTrigger().getMath())).group(1))
        assign = ev.getEventAssignment(0)
        events.append((t_ev, assign.getVariable(), float(re.search(r"\+ ([0-9.eE+-]+)", libsbml.formulaToL3String(assign.getMath())).group(1))))
        breaks.add(t_ev)

    funcs = {"_pw": pw, "exp": math.exp, "sqrt": math.sqrt}

    def rhs(t, y):
        scope = {**consts, **{sid: y[i] for sid, i in pos.items()}, "time": t}
        d = np.zeros(len(ids))
        for code, stoich in reactions:
            try:
                v = eval(code, funcs, scope)
            except ZeroDivisionError:  # 0^-0.3 in a piecewise branch that is not selected at t = 0
                v = 0.0
            for sid, n in stoich.items():
                if sid not in fixed:
                    d[pos[sid]] += n * v
        return d

    t_end = float(max(t_eval))
    cuts = sorted({0.0, t_end} | {b for b in breaks if 0.0 < b < t_end})
    out_t, out_y = [], []
    y = y0.copy()
    for a, b in zip(cuts[:-1], cuts[1:]):
        for t_ev, sid, amount in events:
            if t_ev == a:
                y[pos[sid]] += amount
        grid = np.array([t for t in t_eval if a <= t <= b])
        sol = solve_ivp(rhs, (a, b), y, method="BDF", rtol=rtol, atol=atol, dense_output=True)
        out_t.append(grid)
        out_y.append(sol.sol(grid))
        y = sol.y[:, -1]
    t_all = np.concatenate(out_t)
    y_all = np.concatenate(out_y, axis=1)
    keep = np.concatenate([np.diff(t_all) > 0, [True]])  # a cut time is on two segments: keep the later (post-event) one
    names = {sid: model.getSpecies(sid).getName() for sid in ids}
    return {names[sid]: y_all[pos[sid]][keep] for sid in ids}


def _drug_model_and_bees(competition):
    """Three inhibited reactions (competitive+mixed on S->P, tight-binding on a reversible P->Q,
    competitive on S->W sharing the first enzyme) and a dosed, eliminated drug."""
    from bees.core_edge_model import CoreEdgeModel, RateLawOptions, SpeciesData
    from bees.inhibition import InhibitorSpec

    concs = {"S": 5.0, "P": 0.0, "Q": 0.1, "W": 0.0, "E1": 0.02, "E2": 0.005, "Drug": 0.0, "Cmp": 0.4}
    model = CoreEdgeModel()
    model.rate_law_options = RateLawOptions(enzyme_competition=competition)
    for lab, c in concs.items():
        model.add_core_species(SpeciesData(
            label=lab, concentration=c, initial_concentration=c,
            is_enzyme=lab.startswith("E"), constant=lab.startswith("E") or lab == "Cmp",
        ))

    def rxn(enzyme, reactants, products, km_per, kcat, inhibitors, reversible=False, keq=None):
        r = _make_reaction(enzyme_label=enzyme, reactant_labels=reactants, product_labels=products,
                           stoichiometry={**{x: -1 for x in reactants}, **{x: 1 for x in products}}, kcat=kcat)
        r.kinetics.km = None
        r.kinetics.km_per_substrate = km_per
        r.inhibitors = inhibitors
        r.thermo = SimpleNamespace(source="test", irreversible=not reversible, keq=keq or float("nan"),
                                   kcat_rev=0.0, dgr_prime_kJmol=-5.0, sigma_kJmol=1.0)
        r.template.reversible = reversible
        return r

    drug = [InhibitorSpec("drug", "competitive", 0.08, None, "s"), InhibitorSpec("cmp", "mixed", 0.3, 2.5, "s")]
    for r in (
        rxn("E1", ["S"], ["P"], {"S": 1.5}, 8.0, drug),
        rxn("E1", ["S"], ["W"], {"S": 3.0}, 2.0, [InhibitorSpec("drug", "uncompetitive", 0.2, None, "s")]),
        rxn("E2", ["P"], ["Q"], {"P": 0.7, "Q": 1.2}, 30.0, [InhibitorSpec("drug", "tight_binding", 0.002)],
            reversible=True, keq=40.0),
    ):
        model.add_reaction(r)

    dosing = [
        SimpleNamespace(species="Drug", profile="weibull", dose=0.3, params={"tau": 90.0, "beta": 0.7},
                        t_start=0.0, elimination_k=None, elimination_half_life=150.0),
        SimpleNamespace(species="Drug", profile="bolus", dose=0.05, params=None, t_start=60.0,
                        elimination_k=None, elimination_half_life=None),
    ]
    bees = SimpleNamespace(
        species=[SimpleNamespace(label=k, concentration=v) for k, v in concs.items() if not k.startswith("E")],
        enzymes=[SimpleNamespace(label=k, concentration=v) for k, v in concs.items() if k.startswith("E")],
        dosing=dosing, settings=SimpleNamespace(calibrations=[]),
    )
    return model, bees


@pytest.mark.parametrize("competition", [False, True])
def test_sbml_with_inhibitors_and_dosing_integrates_like_bees(output_dir, competition):
    """The exported SBML, integrated from its own equations, reproduces the BEES trajectory."""
    from bees.dosing import schedules_from_spec
    from bees.simulator import ODESimulator

    pytest.importorskip("libsbml")
    model, bees = _drug_model_and_bees(competition)
    schedules = schedules_from_spec(bees.dosing)
    t_eval = np.arange(0.0, 301.0, 15.0)
    ours = ODESimulator(model, dosing=schedules).simulate(
        end_time=300.0, time_step=15.0, rtol=1e-10, atol=1e-13)
    assert ours.success

    exporter = EnlargerExporter(
        model=model, profiles=[], output_directory=output_dir, logger=MagicMock(),
        reaction_id_by_sig={}, reaction_first_seen_iter={}, reaction_core_enter_iter={},
        reaction_obj_by_sig={}, iteration_summaries=[], bees_object=bees,
    )
    path = exporter.export_sbml(filename="drug.xml")
    sbml = _integrate_sbml(path, t_eval)
    for lab in ("S", "P", "Q", "W", "Drug"):
        assert sbml[lab] == pytest.approx(ours.y[ours.species_labels.index(lab)], rel=2e-4, abs=1e-8), lab
    # the drug really moved the chemistry
    plain_model, _ = _drug_model_and_bees(competition)
    for r in plain_model.core_reactions:
        r.inhibitors = None
    plain = ODESimulator(plain_model).simulate(end_time=300.0, time_step=15.0, rtol=1e-10, atol=1e-13)
    assert ours.y[ours.species_labels.index("Q"), -1] < 0.8 * plain.y[plain.species_labels.index("Q"), -1]


def test_sbml_has_inhibitor_modifiers_parameters_and_dose_reactions(output_dir):
    libsbml = pytest.importorskip("libsbml")
    model, bees = _drug_model_and_bees(False)
    exporter = EnlargerExporter(
        model=model, profiles=[], output_directory=output_dir, logger=MagicMock(),
        reaction_id_by_sig={}, reaction_first_seen_iter={}, reaction_core_enter_iter={},
        reaction_obj_by_sig={}, iteration_summaries=[], bees_object=bees,
    )
    doc = libsbml.readSBML(exporter.export_sbml(filename="drug.xml"))
    m = doc.getModel()
    assert doc.getNumErrors(libsbml.LIBSBML_SEV_ERROR) == 0
    modifiers = {r.getId(): {r.getModifier(k).getSpecies() for k in range(r.getNumModifiers())}
                 for r in (m.getReaction(i) for i in range(m.getNumReactions()))}
    assert {"Drug", "Cmp"} <= modifiers["R1"]
    params = {m.getParameter(i).getId() for i in range(m.getNumParameters())}
    assert any(p.startswith("Ki_R1_Drug") for p in params)
    assert any(p.startswith("alpha_R1_Cmp") for p in params) and any(p.startswith("Ki_R1_Cmp") for p in params)
    assert any(p.startswith("Ki_R3_Drug") for p in params)
    assert m.getReaction("DOSE_1_Drug") is not None and m.getReaction("ELIM_Drug") is not None
    assert m.getReaction("DOSE_2_Drug") is None  # a bolus is not a rate
    assert m.getNumEvents() == 1
    # weibull beta < 1: the first instant (a few 1e-5 of the dose) is part of the initial concentration
    assert m.getSpecies("Drug").getInitialConcentration() == pytest.approx(0.3 * (1 - np.exp(-(1e-3 / 90.0) ** 0.7)))


def test_sbml_bolus_at_time_zero_is_the_initial_concentration(output_dir):
    pytest.importorskip("libsbml")
    model, bees = _drug_model_and_bees(False)
    bees.dosing = [SimpleNamespace(species="Drug", profile="bolus", dose=0.25, params=None, t_start=0.0,
                                   elimination_k=None, elimination_half_life=None)]
    exporter = EnlargerExporter(
        model=model, profiles=[], output_directory=output_dir, logger=MagicMock(),
        reaction_id_by_sig={}, reaction_first_seen_iter={}, reaction_core_enter_iter={},
        reaction_obj_by_sig={}, iteration_summaries=[], bees_object=bees,
    )
    libsbml = pytest.importorskip("libsbml")
    m = libsbml.readSBML(exporter.export_sbml(filename="drug.xml")).getModel()
    assert m.getSpecies("Drug").getInitialConcentration() == pytest.approx(0.25)
    assert m.getNumEvents() == 0 and m.getReaction("ELIM_Drug") is None
