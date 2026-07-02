"""
Small housekeeping helpers that sort files between working directories
based on lookup lists or directory contents, without doing any network or
structural analysis themselves.
"""

from __future__ import annotations

import shutil
from pathlib import Path
from typing import Iterable

from . import config


def sort_known_dimers(pdb_id: str, known_dimer_ids: Iterable[str]) -> None:
    """
    Move files in ``mer/dimer/`` whose stem is a known dimer (from an
    external reference list) into ``mer/dimer/really_dimer/``.

    Equivalent to the original ``is_dimer``'s live code path (the original
    also contained an unreachable ``csv``-based re-derivation of the dimer
    list after an early ``return`` -- that dead branch has been dropped).
    """
    if len(pdb_id) != 4:
        return
    known = set(known_dimer_ids)
    dimer_dir = config.mer_dir(pdb_id) / "dimer"
    if not dimer_dir.is_dir():
        return
    really_dimer_dir = dimer_dir / "really_dimer"

    for entry in dimer_dir.iterdir():
        if entry.stem in known:
            really_dimer_dir.mkdir(parents=True, exist_ok=True)
            shutil.move(str(entry), str(really_dimer_dir))


def promote_orphan_fastas(pdb_id: str) -> None:
    """
    Within ``mer/overlapping/``, promote a candidate's FASTA to
    ``mer/alreadymer/`` if no matching ``<pdb_id>_<candidate>.pdb`` sits
    alongside it (meaning the pairing was never actually modelled).

    Equivalent to the original ``check_fasta_pdb``.
    """
    if len(pdb_id) != 4:
        return
    overlapping_dir = config.mer_dir(pdb_id) / "overlapping"
    if not overlapping_dir.is_dir():
        return

    entries = {p.name for p in overlapping_dir.iterdir()}
    already_mer_dir = config.mer_dir(pdb_id) / "alreadymer"

    for entry in list(overlapping_dir.glob("*.fasta")):
        candidate_id = entry.stem
        expected_pdb = f"{pdb_id}_{candidate_id}.pdb"
        if expected_pdb not in entries:
            already_mer_dir.mkdir(parents=True, exist_ok=True)
            shutil.move(str(entry), str(already_mer_dir))


def reconcile_empty_monomer_dir(pdb_id: str, monomer_root: Path | None = None) -> None:
    """
    If a previously-sorted "confirmed monomer" entry's ``mer/`` directory no
    longer contains any FASTA/PDB evidence, move it back up a level.

    Equivalent to the original ``pdb_not_in_file``.
    """
    if len(pdb_id) != 4:
        return
    monomer_root = monomer_root or (config.bucket_dir("monomer") / "monomer")
    path = monomer_root / pdb_id
    mer_directory = path / "mer"

    if not mer_directory.is_dir():
        shutil.move(str(path), str(config.bucket_dir("monomer")))
        return

    has_evidence = any(e.suffix in (".fasta", ".pdb") for e in mer_directory.iterdir())
    if not has_evidence:
        shutil.move(str(path), str(config.bucket_dir("monomer")))
