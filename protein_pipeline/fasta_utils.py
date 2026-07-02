"""
FASTA retrieval and coarse alignment-quality filtering.

Replaces ``get_fasta``, ``alignedResidues`` and ``not_monomer`` from the
original script.
"""

from __future__ import annotations

import logging
import shutil
from pathlib import Path

from . import config, http_client, vast

logger = logging.getLogger(__name__)

_SKIP_NAMES = {"no_aligne", "monomer", "false_mer"}


def fetch_query_fasta(pdb_id: str) -> Path | None:
    """
    Ensure ``<pdb_id>/<pdb_id>.fasta`` exists, fetching it from RCSB if not.

    Mirrors ``get_fasta``, which downloaded a FASTA for every entry sitting
    in a PDB's working directory (the query plus any not-yet-fetched
    neighbours) rather than just the query; that broader behaviour is
    exposed separately as :func:`fetch_missing_fastas`.
    """
    if pdb_id in ("oligState", "no_homomer", "no"):
        return None
    target = config.pdb_dir(pdb_id) / f"{pdb_id}.fasta"
    if target.exists():
        return target
    return http_client.download_to_file(http_client.fetch_rcsb_fasta, pdb_id, target)


def fetch_missing_fastas(pdb_id: str) -> None:
    """Fetch a FASTA for every file currently in ``<pdb_id>/`` lacking one."""
    if pdb_id in ("oligState", "no_homomer", "no"):
        return
    path = config.pdb_dir(pdb_id)
    for entry in path.iterdir():
        stem = entry.name.split("_")[0]
        http_client.download_to_file(http_client.fetch_rcsb_fasta, stem, path / f"{stem}.fasta")


def _read_sequence(fasta_path: Path) -> str:
    """Return the sequence line (line 2) of a FASTA file, or ``''`` on failure."""
    try:
        return fasta_path.read_text().split("\n")[1]
    except (FileNotFoundError, IndexError):
        return ""


def filter_aligned_neighbours(pdb_id: str) -> None:
    """
    Discard candidate homomer partners whose reported aligned-residue ratio
    (relative to the query length) falls outside
    :data:`config.ALIGNED_RESIDUE_RATIO_RANGE`, or whose sequence length
    ratio to the query is likewise out of range.

    Equivalent to the original ``alignedResidues``. Rejected candidates are
    moved into ``<pdb_id>/mer/no_aligne/``.
    """
    if len(pdb_id) != 4:
        return

    path = config.pdb_dir(pdb_id)
    query_sequence = _read_sequence(path / f"{pdb_id}.fasta")
    if not query_sequence:
        logger.warning("No query sequence for %s; skipping alignment filter", pdb_id)
        return

    no_aligne_dir = config.mer_dir(pdb_id) / "no_aligne"
    no_aligne_dir.mkdir(parents=True, exist_ok=True)

    neighbours = {n.pdb_id: n.aligned_residues for n in vast.load_neighbours(pdb_id)}
    lo, hi = config.ALIGNED_RESIDUE_RATIO_RANGE

    for entry in list(config.mer_dir(pdb_id).iterdir()):
        if entry.name in _SKIP_NAMES or entry.suffix == ".xml" or not entry.name.endswith(".fasta"):
            continue

        neighbour_id = entry.stem
        candidate_sequence = _read_sequence(entry)

        if not candidate_sequence:
            _reject(pdb_id, neighbour_id, no_aligne_dir)
            continue

        aligned_ratio = neighbours.get(neighbour_id, 0) / len(query_sequence)
        length_ratio = len(candidate_sequence) / len(query_sequence)
        if not (lo <= aligned_ratio <= hi) or not (lo <= length_ratio <= hi):
            _reject(pdb_id, neighbour_id, no_aligne_dir)

    if len(list(config.mer_dir(pdb_id).iterdir())) == 1:
        shutil.move(str(path), str(config.bucket_dir("no_aligne")))


def _reject(pdb_id: str, neighbour_id: str, no_aligne_dir: Path) -> None:
    """Move a rejected candidate's fasta into ``no_aligne`` and drop its xml."""
    fasta_path = config.mer_dir(pdb_id) / f"{neighbour_id}.fasta"
    xml_path = config.mer_dir(pdb_id) / f"{neighbour_id}.xml"
    if fasta_path.exists():
        shutil.move(str(fasta_path), str(no_aligne_dir))
    if xml_path.exists():
        xml_path.unlink()


def split_out_true_monomers(pdb_id: str) -> None:
    """
    Move candidate neighbours that VAST+ itself reports as single-molecule
    (``numberOfMolecules == 1``) into ``<pdb_id>/mer/monomer/``.

    Equivalent to the original ``not_monomer`` (which, despite the name,
    filters out entries that *are* monomeric).
    """
    if pdb_id in ("oligState", "no_homomer", "no", "no_aligne", "monomer"):
        return

    monomer_dir = config.mer_dir(pdb_id) / "monomer"
    monomer_dir.mkdir(parents=True, exist_ok=True)

    molecule_counts = {n.pdb_id: n.number_of_molecules for n in vast.load_neighbours(pdb_id)}

    for entry in list(config.mer_dir(pdb_id).iterdir()):
        if entry.name in ("no_aligne", "monomer"):
            continue
        neighbour_id = entry.name.split(".")[0]
        if molecule_counts.get(neighbour_id) == 1:
            shutil.move(str(entry), str(monomer_dir))

    if len(list(config.mer_dir(pdb_id).iterdir())) == 2:
        shutil.move(str(config.pdb_dir(pdb_id)), str(config.bucket_dir("monomer")))
