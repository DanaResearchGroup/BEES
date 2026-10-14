"""Shared helpers for the BEES rule layer (laws + calibrations)."""

from __future__ import annotations

from typing import Optional

from bees.cofactors import COFACTORS_ALWAYS_AVAILABLE

_CO2_SMILES = "O=C=O"
_CO2_LABELS = frozenset({"carbon dioxide", "co2", "carbon-dioxide", "co₂"})

def detect_acyl_chain_length(label: str, smiles: Optional[str] = None) -> Optional[int]:
    """Acyl carbon count (incl. carbonyl) from structure, or None.

    The label only excludes buffered cofactors and CO2; no SMILES means no chain length.
    """
    lab = (label or "").lower().strip()
    if lab in COFACTORS_ALWAYS_AVAILABLE or lab in _CO2_LABELS:
        return None
    if not smiles:
        return None
    return _acyl_chain_from_smiles(smiles)

def _acyl_anchor(mol) -> tuple[Optional[int], Optional[int]]:
    """(carbonyl C, S/O atom it is bound to): a thioester if any, else an ester/acid.

    Thioester first so malonyl thioesters anchor on the carrier, not on the free carboxyl.
    """
    for want in ("S", "O"):
        for atom in mol.GetAtoms():
            if atom.GetSymbol() != "C":
                continue
            has_carbonyl = False
            hetero = None
            for nb in atom.GetNeighbors():
                order = mol.GetBondBetweenAtoms(atom.GetIdx(), nb.GetIdx()).GetBondTypeAsDouble()
                if nb.GetSymbol() == "O" and order == 2.0:
                    has_carbonyl = True
                elif nb.GetSymbol() == want and order == 1.0:
                    hetero = nb.GetIdx()
            if has_carbonyl and hetero is not None:
                return atom.GetIdx(), hetero
    return None, None

def _acyl_chain_from_smiles(smiles: str) -> Optional[int]:
    """Count acyl carbons from SMILES (thioester/carboxyl anchor)."""
    try:
        from rdkit import Chem
    except Exception:
        return None
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None

    anchor_idx, _ = _acyl_anchor(mol)
    if anchor_idx is None:
        return None

    seen: set[int] = set()
    stack = [anchor_idx]
    while stack:
        idx = stack.pop()
        if idx in seen:
            continue
        seen.add(idx)
        atom = mol.GetAtomWithIdx(idx)
        for nb in atom.GetNeighbors():
            if nb.GetSymbol() == "C" and nb.GetIdx() not in seen:
                stack.append(nb.GetIdx())
    return len(seen) if seen else None

def detect_acyl_carrier(smiles: Optional[str]) -> Optional[str]:
    """Canonical SMILES of what the acyl carbonyl is bound to (e.g. 'O' for a free acid), or None.

    Charge and stereo are dropped so one carrier always gives one key.
    """
    if not smiles:
        return None
    try:
        from rdkit import Chem
        from rdkit.Chem.MolStandardize import rdMolStandardize
    except Exception:
        return None
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None
    mol = rdMolStandardize.Uncharger().uncharge(mol)
    anchor_idx, hetero_idx = _acyl_anchor(mol)
    if anchor_idx is None:
        return None
    seen: set[int] = set()
    stack = [hetero_idx]
    while stack:
        idx = stack.pop()
        if idx in seen:
            continue
        seen.add(idx)
        for nb in mol.GetAtomWithIdx(idx).GetNeighbors():
            if nb.GetIdx() != anchor_idx and nb.GetIdx() not in seen:
                stack.append(nb.GetIdx())
    return Chem.MolFragmentToSmiles(mol, atomsToUse=sorted(seen), isomericSmiles=False)

def _ec_norm(ec) -> Optional[str]:
    """'EC 3.1.2.14' -> '3.1.2.14'."""
    if not ec:
        return None
    s = str(ec).strip().upper()
    if s.startswith("EC "):
        s = s[3:].strip()
    return s or None

def _reaction_acyl_substrate_n(reaction) -> Optional[int]:
    """Longest acyl-chain length among substrates, or None."""
    kin = getattr(reaction, "kinetics", None)
    stoich = getattr(reaction, "stoichiometry", None) or {}
    smiles_map = (getattr(kin, "compound_smiles", None) or {}) if kin is not None else {}
    best: Optional[int] = None
    for lab, coeff in stoich.items():
        if coeff >= 0:
            continue
        if str(lab).lower().strip() in COFACTORS_ALWAYS_AVAILABLE:
            continue
        n = detect_acyl_chain_length(lab, smiles_map.get(lab))
        if n is not None and (best is None or n > best):
            best = n
    return best
