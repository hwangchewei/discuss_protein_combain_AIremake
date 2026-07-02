"""
Central configuration for the protein_pipeline package.

The original script scattered hard-coded Windows paths
(``D:\\programming_language\\python\\monomer_json``, ``E:\\monomer_json`` ...)
throughout every function. Here all paths are derived from a single
``PROJECT_ROOT`` that can be overridden with the ``PROTEIN_PIPELINE_ROOT``
environment variable, and every path is a ``pathlib.Path`` so the code works
the same way on Windows, macOS, and Linux.
"""

from __future__ import annotations

import os
from pathlib import Path

# ---------------------------------------------------------------------------
# Root locations
# ---------------------------------------------------------------------------

#: Root directory for all pipeline data. Override with the
#: PROTEIN_PIPELINE_ROOT environment variable.
PROJECT_ROOT = Path(os.environ.get("PROTEIN_PIPELINE_ROOT", "./monomer_json")).resolve()

#: Where per-PDB working directories (``<root>/<pdb_id>/``) live.
DATA_DIR = PROJECT_ROOT

#: Bucket directories the pipeline sorts finished/rejected entries into.
BUCKET_NAMES = (
    "no_homomer",
    "oligState",
    "no_aligne",
    "false_mer",
    "alreadymer",
    "monomer",
    "overlapping",
    "interface_error",
)


def bucket_dir(name: str) -> Path:
    """Return (and ensure the existence of) one of the sorting-bucket dirs."""
    path = DATA_DIR / name
    path.mkdir(parents=True, exist_ok=True)
    return path


def pdb_dir(pdb_id: str) -> Path:
    """Working directory for a single PDB entry, e.g. ``<root>/1ABC``."""
    return DATA_DIR / pdb_id


def mer_dir(pdb_id: str) -> Path:
    """Directory holding candidate homomer partners for a given entry."""
    return pdb_dir(pdb_id) / "mer"


# ---------------------------------------------------------------------------
# External tools
# ---------------------------------------------------------------------------

#: Path (or bare command name, if it's on PATH) to the chromedriver binary.
CHROMEDRIVER_PATH = os.environ.get("CHROMEDRIVER_PATH", "chromedriver")

#: Path (or bare command name) to the clustalw2 executable.
CLUSTALW_EXE = os.environ.get("CLUSTALW_EXE", "clustalw2")

# ---------------------------------------------------------------------------
# Network
# ---------------------------------------------------------------------------

USER_AGENT_STRING = (
    "Mozilla/5.0 (Windows NT 10.0; Win64; x64) AppleWebKit/537.36 "
    "(KHTML, like Gecko) Chrome/120.0.0.0 Safari/537.36"
)
REQUEST_HEADERS = {"User-Agent": USER_AGENT_STRING}
REQUEST_TIMEOUT = 30  # seconds

# ---------------------------------------------------------------------------
# Filtering thresholds (previously magic numbers scattered in the script)
# ---------------------------------------------------------------------------

#: Maximum RMSD (Angstrom) for a VAST+ neighbour to be considered a match.
MAX_NEIGHBOUR_RMSD = 1.7

#: Minimum sequence identity (0-1) for a VAST+ neighbour to be considered a match.
MIN_NEIGHBOUR_SEQUENCE_IDENTITY = 0.2

#: Aligned-residue-count ratio window (relative to the query length) used to
#: decide whether an alignment is acceptable.
ALIGNED_RESIDUE_RATIO_RANGE = (0.9, 1.1)

#: Interfacial solvation energy (kcal/mol) threshold used to flag a genuine
#: protein-protein interface in PISA data.
MAX_INTERFACE_SOLVATION_ENERGY = -8.0
