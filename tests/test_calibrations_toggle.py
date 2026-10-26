"""Unit tests for RuleRegistry.configure_calibrations and per-run registry state.

Verifies that settings.calibrations enables exactly the named calibrations,
never a physics law or an unknown name, that state does not leak across
calls (the global-singleton reset), and that calibration observables follow
the enabled set.
"""
from __future__ import annotations

import pytest

from bees.rules import Reference, Rule, RuleRegistry


def _rule(name: str, kind: str, enabled: bool = True) -> Rule:
    return Rule(
        name=name,
        description="test",
        reference=Reference(authors=("Test",), title="t", year="2024"),
        reference_type="theoretical",
        kind=kind,
        enabled=enabled,
    )


def _registry() -> RuleRegistry:
    reg = RuleRegistry()
    reg.register(_rule("a_law", "law", enabled=True))
    reg.register(_rule("cal_one", "calibration", enabled=False))
    reg.register(_rule("cal_two", "calibration", enabled=False))
    return reg


def _enabled(reg):
    return {r.name: r.enabled for r in reg}


def test_enables_exactly_named_calibration():
    reg = _registry()
    reg.configure_calibrations(["cal_one"])
    assert _enabled(reg) == {"a_law": True, "cal_one": True, "cal_two": False}


def test_empty_turns_previous_calibrations_off_law_untouched():
    reg = _registry()
    reg.configure_calibrations(["cal_one", "cal_two"])
    # A subsequent run that lists nothing must turn them back off (no leak across calls).
    reg.configure_calibrations([])
    assert _enabled(reg) == {"a_law": True, "cal_one": False, "cal_two": False}


def test_none_is_treated_as_empty():
    reg = _registry()
    reg.configure_calibrations(["cal_one", "cal_two"])
    reg.configure_calibrations(None)
    assert _enabled(reg) == {"a_law": True, "cal_one": False, "cal_two": False}


def test_law_name_rejected():
    reg = _registry()
    with pytest.raises(ValueError, match="always-on physics law"):
        reg.configure_calibrations(["a_law"])


def test_unknown_name_rejected_lists_available():
    reg = _registry()
    with pytest.raises(ValueError, match=r"Unknown calibration 'nope'.*Available calibrations: .*cal_one.*cal_two"):
        reg.configure_calibrations(["nope"])


def test_reset_baselines_forgets_a_previous_models_values():
    from types import SimpleNamespace

    class Double(Rule):
        def apply(self, reaction):
            reaction.kinetics.kcat *= 2

    def rxn(kcat):
        return SimpleNamespace(
            enzyme_label="E", reactant_labels=["S"], product_labels=["P"],
            kinetics=SimpleNamespace(kcat=kcat, km=1.0, km_per_substrate=None), thermo=None,
        )

    reg = RuleRegistry()
    reg.register(Double(name="double", description="t",
                        reference=Reference(authors=("T",), title="t", year="2024"),
                        reference_type="theoretical"))
    first = rxn(10.0)
    reg.apply_all([first])
    assert first.kinetics.kcat == 20.0

    reg.reset_baselines()
    second = rxn(50.0)  # same signature, new model
    reg.apply_all([second])
    assert second.kinetics.kcat == 100.0


def test_observable_emitted_only_when_its_calibration_is_enabled():
    reg = _registry()
    reg.register_observable(lambda ids: ("obs", "Obs", " + ".join(ids.values())),
                            calibrations=["cal_one"])
    ids = {"A": "A_id", "B": "B_id"}
    assert reg.observables([], ids) == []
    assert reg.observables(["cal_two"], ids) == []
    assert reg.observables(["cal_one"], ids) == [("obs", "Obs", "A_id + B_id")]


def test_validation_failure_leaves_state_unchanged():
    # An invalid request should raise before flipping any enabled flags.
    reg = _registry()
    reg.configure_calibrations(["cal_one"])
    with pytest.raises(ValueError):
        reg.configure_calibrations(["cal_two", "a_law"])
    # cal_one still enabled, cal_two still disabled (no partial application).
    assert reg.by_name("cal_one").enabled
    assert not reg.by_name("cal_two").enabled
