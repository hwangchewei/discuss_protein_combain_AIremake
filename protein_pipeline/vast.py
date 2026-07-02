"""
NCBI VAST+ structural-neighbour discovery.

VAST+ compares a query structure against the PDB and returns a JSON payload
listing structural neighbours together with RMSD, sequence identity, and
oligomeric-state metadata. This module downloads that payload and filters
the neighbour list down to plausible homomer partners.

Replaces ``json_download``, ``no_selection``, ``monomer_fasta`` and
``oligState`` from the original script. The original ``json_download`` drove
a full Chrome browser just to click a "no selection" confirmation button
before the JSON became downloadable; that interactive step is preserved in
:func:`download_vastplus_json_via_browser` for the (site-dependent) cases
where the plain HTTP endpoint isn't sufficient, but :func:`fetch_vastplus_json`
should be tried first since it avoids the Selenium dependency entirely.
"""

from __future__ import annotations

import json
import logging
from dataclasses import dataclass
from pathlib import Path
from time import sleep
from typing import List

from selenium.webdriver.common.by import By

from . import browser, config, http_client

logger = logging.getLogger(__name__)


@dataclass
class Neighbour:
    """A single VAST+ structural neighbour entry."""

    pdb_id: str
    rmsd: float
    sequence_identity: float
    aligned_residues: int
    number_of_molecules: int

    @classmethod
    def from_json(cls, entry: dict) -> "Neighbour":
        return cls(
            pdb_id=entry["pdbId"],
            rmsd=float(entry["rmsd"]),
            sequence_identity=float(entry["sequenceIdentity"]),
            aligned_residues=int(entry["alignedResidues"]),
            number_of_molecules=int(entry["numberOfMolecules"]),
        )


def download_vastplus_json_via_browser(pdb_id: str, download_dir: Path) -> Path | None:
    """
    Drive a headed Chrome session to trigger the VAST+ JSON download for
    ``pdb_id``, for cases where the plain HTTP endpoint requires the
    interactive "no selection" confirmation click first.

    Returns the path the browser saved the file to, or ``None`` if the PDB
    id had no VAST+ data (recorded to ``no_data.txt`` as in the original).
    """
    with browser.chrome_driver(headless=False, download_dir=str(download_dir)) as bro:
        bro.get(f"https://www.ncbi.nlm.nih.gov/Structure/vastplus/vastplus.cgi?uid={pdb_id}")
        try:
            bro.find_element(By.XPATH, "/html/body/div[1]/div/div[3]/div[1]/div[2]/button[2]").click()
        except Exception:
            (download_dir / "no_data.txt").open("a").write(pdb_id + "\n")
            return None
        sleep(2)
    return download_dir / f"{pdb_id}_vastplus.json"


def download_vastplus_json_via_browser_or_http(pdb_id: str) -> Path | None:
    """
    Fetch and save the VAST+ JSON for ``pdb_id``, trying the plain HTTP
    endpoint first and only falling back to the interactive browser flow
    (:func:`download_vastplus_json_via_browser`) if that fails.

    Suitable for use as the target function of a ``multiprocessing.Pool``,
    mirroring how the original script's ``main()`` parallelised
    ``json_download`` across a pool of worker processes.
    """
    destination_dir = config.pdb_dir(pdb_id)
    destination = destination_dir / f"{pdb_id}_vastplus.json"
    try:
        http_client.download_to_file(http_client.fetch_vastplus_json, pdb_id, destination)
        return destination
    except Exception:
        logger.info("Plain HTTP VAST+ fetch failed for %s, falling back to browser", pdb_id)
        return download_vastplus_json_via_browser(pdb_id, destination_dir)


def load_neighbours(pdb_id: str) -> List[Neighbour]:
    """Load and parse the ``<pdb_id>_vastplus.json`` file for an entry."""
    json_path = config.pdb_dir(pdb_id) / f"{pdb_id}_vastplus.json"
    payload = json.loads(json_path.read_text())
    return [Neighbour.from_json(entry) for entry in payload.get("neighbors", [])]


def query_oligomeric_state(pdb_id: str) -> str | None:
    """Return the VAST+-reported oligomeric state of the query structure, if present."""
    json_path = config.pdb_dir(pdb_id) / f"{pdb_id}_vastplus.json"
    payload = json.loads(json_path.read_text())
    return payload.get("query", {}).get("oligState")


def is_plausible_homomer(neighbour: Neighbour, monomer_ids: set, homomer_ids: set) -> bool:
    """
    Decide whether a VAST+ neighbour looks like a genuine homomer partner.

    Mirrors the filter in the original ``monomer_fasta``: the neighbour must
    not itself be a known monomer, must appear in the known-homomer list,
    and must pass the RMSD / sequence-identity thresholds.
    """
    pdb_id = neighbour.pdb_id.lower()
    return (
        pdb_id not in monomer_ids
        and pdb_id in homomer_ids
        and neighbour.rmsd < config.MAX_NEIGHBOUR_RMSD
        and neighbour.sequence_identity > config.MIN_NEIGHBOUR_SEQUENCE_IDENTITY
    )


def is_candidate_neighbour(neighbour: Neighbour, monomer_ids: set) -> bool:
    """Looser filter used to decide whether a neighbour's own VAST+ JSON should be kept."""
    return (
        neighbour.pdb_id.lower() not in monomer_ids
        and neighbour.rmsd < config.MAX_NEIGHBOUR_RMSD
        and neighbour.sequence_identity > config.MIN_NEIGHBOUR_SEQUENCE_IDENTITY
    )


def _read_id_list(path: Path) -> set:
    """Parse a whitespace-delimited id list file (col 2 = PDB id), lower-cased."""
    lines = path.read_text().splitlines()
    return {line.split()[1].lower() for line in lines if len(line.split()) > 1}


def collect_candidate_fastas(
    pdb_id: str,
    monomer_list_path: Path,
    homomer_list_path: Path,
    *,
    require_homomer: bool = True,
) -> int:
    """
    Given a downloaded VAST+ JSON for ``pdb_id``, fetch FASTA files for
    plausible homomer-partner neighbours into ``<pdb_id>/mer/``.

    Combines ``no_selection`` (``require_homomer=False``, fetch everything)
    and ``monomer_fasta`` (``require_homomer=True``, apply the RMSD/identity/
    known-homomer filter) from the original script, which were otherwise
    near-identical copies of each other.

    :returns: number of FASTA files written.
    """
    monomer_ids = _read_id_list(monomer_list_path)
    homomer_ids = _read_id_list(homomer_list_path) if require_homomer else set()

    mer_directory = config.mer_dir(pdb_id)
    mer_directory.mkdir(parents=True, exist_ok=True)

    written = 0
    for neighbour in load_neighbours(pdb_id):
        keep = (
            is_plausible_homomer(neighbour, monomer_ids, homomer_ids)
            if require_homomer
            else True
        )
        if keep:
            http_client.download_to_file(
                http_client.fetch_rcsb_fasta,
                neighbour.pdb_id,
                mer_directory / f"{neighbour.pdb_id}.fasta",
            )
            written += 1

        if require_homomer and is_candidate_neighbour(neighbour, monomer_ids):
            candidate_json = config.pdb_dir(pdb_id) / f"{neighbour.pdb_id}_vastplus.json"
            if candidate_json.exists():
                candidate_json.rename(config.pdb_dir(pdb_id) / candidate_json.name)

    if written == 0:
        import shutil
        shutil.move(str(config.pdb_dir(pdb_id)), str(config.bucket_dir("no_homomer")))
    return written


def sort_by_oligomeric_state(pdb_id: str) -> None:
    """Move an entry into the ``oligState`` bucket if it isn't reported as monomeric."""
    if pdb_id in ("oligState", "no_homomer", "no"):
        return
    state = query_oligomeric_state(pdb_id)
    if state != "monomeric":
        import shutil
        shutil.move(str(config.pdb_dir(pdb_id)), str(config.bucket_dir("oligState")))
