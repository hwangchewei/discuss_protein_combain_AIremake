"""
Amino-acid codes and physicochemical classification.

Consolidates the several near-duplicate copies of ``aa_codes`` /
``Pcharge`` / ``Ncharge`` / ``Hydrophilic`` / ``Hydrophobic`` that were
defined independently in both ``discuss_protein_combain.py`` and
``__init__.py`` (the ``biopython`` helper module), and replaces the ~40-way
``if/elif`` chains used to tally single-letter amino-acid counts with a
simple ``collections.Counter``-based helper.
"""

from __future__ import annotations

from collections import Counter
from enum import Enum
from typing import Iterable

#: Three-letter -> one-letter amino acid code lookup.
THREE_TO_ONE = {
    "ALA": "A", "CYS": "C", "ASP": "D", "GLU": "E",
    "PHE": "F", "GLY": "G", "HIS": "H", "LYS": "K",
    "ILE": "I", "LEU": "L", "MET": "M", "ASN": "N",
    "PRO": "P", "GLN": "Q", "ARG": "R", "SER": "S",
    "THR": "T", "VAL": "V", "TYR": "Y", "TRP": "W",
}

#: One-letter -> three-letter amino acid code lookup.
ONE_TO_THREE = {one: three for three, one in THREE_TO_ONE.items()}

# Backwards-compatible aliases matching the original script's names.
aa_codes = THREE_TO_ONE
aa_code_reverse = ONE_TO_THREE


class ResidueClass(str, Enum):
    """Coarse physicochemical classification of an amino acid."""

    POSITIVE = "P"
    NEGATIVE = "N"
    HYDROPHILIC = "Lic"
    HYDROPHOBIC = "Bic"
    OTHER = "other"


#: Positively charged (basic) residues.
POSITIVE_CHARGE = {"HIS", "LYS", "ARG"}
#: Negatively charged (acidic) residues.
NEGATIVE_CHARGE = {"ASP", "GLU"}
#: Polar/uncharged residues.
HYDROPHILIC = {"CYS", "ASN", "GLN", "SER", "THR", "TYR"}
#: Non-polar residues.
HYDROPHOBIC = {"ALA", "PHE", "ILE", "LEU", "MET", "PRO", "VAL", "TRP", "GLY"}

# Backwards-compatible aliases matching the original script's names.
Pcharge = sorted(POSITIVE_CHARGE)
Ncharge = sorted(NEGATIVE_CHARGE)
Hydrophilic = sorted(HYDROPHILIC)
Hydrophobic = sorted(HYDROPHOBIC)


def classify_residue(three_letter_code: str) -> ResidueClass:
    """Classify a three-letter residue code into a :class:`ResidueClass`."""
    if three_letter_code in POSITIVE_CHARGE:
        return ResidueClass.POSITIVE
    if three_letter_code in NEGATIVE_CHARGE:
        return ResidueClass.NEGATIVE
    if three_letter_code in HYDROPHILIC:
        return ResidueClass.HYDROPHILIC
    if three_letter_code in HYDROPHOBIC:
        return ResidueClass.HYDROPHOBIC
    return ResidueClass.OTHER


def composition(residues: Iterable[str], *, one_letter: bool = False) -> Counter:
    """
    Tally residue classes for a sequence of residue codes.

    Replaces the repeated 20-branch ``if/elif`` chains in the original
    script (see ``interface_amino_change`` / ``protein_type``) with a single
    ``Counter`` over :class:`ResidueClass` values.

    :param residues: an iterable of one- or three-letter residue codes.
    :param one_letter: set True if ``residues`` uses one-letter codes.
    """
    counts: Counter = Counter()
    for res in residues:
        code = ONE_TO_THREE.get(res, res) if one_letter else res
        counts[classify_residue(code)] += 1
    return counts


def letter_counts(residues: Iterable[str]) -> Counter:
    """Tally raw one-letter amino-acid frequencies (A, C, D, E, ...)."""
    return Counter(r for r in residues if r in ONE_TO_THREE)


#: Empirical per-residue weighting used by ``structure_surface.SurfaceAnalyzer``
#: to score candidate oligomerisation surface patches. Kept as data (rather
#: than duplicated in two files as in the original script).
SURFACE_WEIGHT_FUNCTION = {
    "CYS": 1.12, "ASN": -0.05, "GLN": 0.16, "SER": -0.39, "THR": 0.46,
    "TYR": 0.01, "ASP": -0.28, "GLU": -0.30, "HIS": 0.20, "LYS": 0.07,
    "ARG": 0.02, "ALA": -0.06, "PHE": 0.43, "ILE": 0.70, "LEU": 1.01,
    "MET": 0.86, "PRO": 0.83, "VAL": 0.38, "TRP": 1.69, "GLY": -0.25,
}

#: A second, refitted version of :data:`SURFACE_WEIGHT_FUNCTION`; this is the
#: variant actually used by the live (non-commented-out) scoring path in the
#: original ``biopython.surface_area_amount``.
SURFACE_WEIGHT_FUNCTION_2 = {
    "CYS": 1.08, "ASN": -0.16, "GLN": 0.10, "SER": -0.45, "THR": 0.35,
    "TYR": -0.03, "ASP": -0.34, "GLU": -0.37, "HIS": 0.05, "LYS": 0.00,
    "ARG": -0.04, "ALA": -0.17, "PHE": 0.30, "ILE": 0.49, "LEU": 0.85,
    "MET": 0.66, "PRO": 0.59, "VAL": 0.26, "TRP": 1.46, "GLY": -0.43,
}

#: Charged residues (union of positive + negative) as used for surface
#: composition scoring.
CHARGED = POSITIVE_CHARGE | NEGATIVE_CHARGE
