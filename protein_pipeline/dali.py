"""
DALI server automation.

The DALI server (ekhidna2.biocenter.helsinki.fi) superposes a query chain
onto every chain of a target PDB entry and reports per-chain coordinates.
The original script (``dali`` / ``dali_start``) used this to reconstruct a
"pseudo-assembly" PDB file for a candidate homomer partner: it uploads the
query PDB, forces it to align to chain A of the target, then walks through
each aligned chain, extracting the transformed coordinates from a text
result page and concatenating them (with each chain's atom records relabeled
by chain letter) into a single PDB file.

This is by far the most fragile, page-structure-dependent part of the
pipeline; kept close to the original logic (including its retry-by-recursion
pattern, bounded to 3 attempts) but with duplicated file paths and the
raw-string-formatted paths replaced by ``pathlib``.
"""

from __future__ import annotations

import logging
import shutil
from pathlib import Path
from time import sleep

from selenium.webdriver.common.by import By

from . import browser, config, http_client

logger = logging.getLogger(__name__)

DALI_URL = "http://ekhidna2.biocenter.helsinki.fi/dali/"
MAX_ATTEMPTS = 3
POLL_INTERVAL_SECONDS = 5
MAX_POLL_ATTEMPTS = 120


def process_all_candidates(pdb_id: str) -> None:
    """Run :func:`superpose_candidate` for every FASTA candidate in ``mer/``."""
    if len(pdb_id) != 4:
        return
    mer_directory = config.mer_dir(pdb_id)
    for fasta_path in mer_directory.glob("*.fasta"):
        superpose_candidate(pdb_id, fasta_path.stem)

    if not any(entry.suffix for entry in mer_directory.iterdir()):
        shutil.move(str(config.pdb_dir(pdb_id)), str(config.bucket_dir("monomer") / "monomer"))


def superpose_candidate(pdb_id: str, candidate_id: str, attempt: int = 1) -> None:
    """
    Superpose ``pdb_id`` chain A onto every chain of ``candidate_id`` via the
    DALI server, writing a reconstructed assembly PDB file. Equivalent to
    the original ``dali_start``.
    """
    mer_directory = config.mer_dir(pdb_id)

    if attempt > MAX_ATTEMPTS:
        logger.error("DALI superposition failed after %d attempts for %s/%s", MAX_ATTEMPTS, pdb_id, candidate_id)
        (mer_directory / f"{candidate_id}.txt").write_text(f"error{attempt}")
        return

    candidate_pdb_path = mer_directory / f"{candidate_id}.pdb"
    http_client.download_to_file(http_client.fetch_pdb_file, candidate_id, candidate_pdb_path)
    crystal_line = next(
        (line for line in candidate_pdb_path.read_text().split("\n") if "CRYST1" in line), ""
    )

    assembly_path = mer_directory / f"{pdb_id}_{candidate_id}.pdb"
    try:
        chain_count = _run_dali_session(pdb_id, candidate_pdb_path, assembly_path, crystal_line)
    except _DaliRetry:
        return superpose_candidate(pdb_id, candidate_id, attempt + 1)

    if chain_count == 1:
        # Only chain A matched: candidate is itself monomeric.
        monomer_dir = mer_directory / "monomer"
        monomer_dir.mkdir(parents=True, exist_ok=True)
        for suffix, name in (
            (".xml", candidate_id), (".fasta", candidate_id),
            (".pdb", f"{pdb_id}_{candidate_id}"), (".pdb", candidate_id),
        ):
            path = mer_directory / f"{name}{suffix}"
            if path.exists():
                shutil.move(str(path), str(monomer_dir))


class _DaliRetry(Exception):
    """Raised internally to signal the DALI session should be retried."""


def _run_dali_session(pdb_id: str, candidate_pdb_path: Path, assembly_path: Path, crystal_line: str) -> int:
    """Drive the DALI web form; returns the number of aligned chains found."""
    with browser.chrome_driver(headless=True) as bro:
        bro.get(DALI_URL)
        bro.find_element(By.XPATH, '//*[@id="tabs"]/ul/li[5]').click()
        bro.find_element(By.XPATH, '//*[@id="tabs-2"]/div/form/input[2]').send_keys(str(candidate_pdb_path))
        bro.find_element(By.XPATH, '//*[@id="tabs-2"]/div/form/div[4]/div[1]/input[1]').click()
        bro.find_element(By.XPATH, '//*[@id="tabs-2"]/div/form/div[4]/div[2]/input[1]').send_keys(f"{pdb_id}A")
        bro.find_element(By.XPATH, '//*[@id="tabs-2"]/div/form/div[7]/input[2]').click()

        chain_rows = _wait_for_results(bro)

        with assembly_path.open("w") as fp:
            fp.write(crystal_line + "\n")
            if len(chain_rows) == 1:
                fp.write("END\n")
                return 1
            for chain_index in range(len(chain_rows)):
                _extract_chain(bro, fp, chain_index)
            fp.write("END\n")
        return len(chain_rows)


def _wait_for_results(bro):
    for _ in range(MAX_POLL_ATTEMPTS):
        rows = bro.find_elements(By.XPATH, "/html/body/ul/table/tbody/tr")
        if rows:
            return rows
        sleep(POLL_INTERVAL_SECONDS)
    raise _DaliRetry("Timed out waiting for DALI results")


def _extract_chain(bro, fp, chain_index: int) -> None:
    """Open one aligned-chain result, extract its ATOM records, and append them to ``fp``."""
    try:
        bro.find_element(By.XPATH, f"/html/body/ul/table/tbody/tr[{chain_index + 1}]/td[1]/a").click()
        bro.find_element(By.XPATH, "/html/body/form/pre/a[2]").click()
    except Exception:
        raise _DaliRetry(f"Could not open DALI chain result #{chain_index}")

    text_lines = bro.find_element(By.XPATH, "/html/body/pre").text.split("\n")
    chain_letter = chr(65 + chain_index)

    for line in text_lines:
        if " A " in line and "ATOM" in line:
            if chain_letter != "A":
                line = line.replace(" A ", f" {chain_letter} ")
            fp.write(line + "\n")
    fp.write("TER\n")

    sleep(0.3)
    bro.back()
    sleep(0.3)
    bro.back()
