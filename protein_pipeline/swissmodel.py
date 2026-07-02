"""
SWISS-MODEL homology modelling automation.

For each candidate homomer partner, build a homology model of the query
sequence threaded onto the candidate's structure via the SWISS-MODEL
interactive web service. Equivalent to the original ``swiss_model``.
"""

from __future__ import annotations

import logging
from pathlib import Path
from time import sleep

from selenium.webdriver.common.by import By
from selenium.webdriver.support import expected_conditions as EC
from selenium.webdriver.support.ui import WebDriverWait

from . import alignment, browser, config

logger = logging.getLogger(__name__)

SWISSMODEL_URL = "https://swissmodel.expasy.org/interactive#alignment"
INITIAL_WAIT_SECONDS = 100
POLL_INTERVAL_SECONDS = 60


def build_models_for_all_candidates(pdb_id: str) -> None:
    """Build a SWISS-MODEL homology model against every ``mer/*.fasta`` candidate."""
    path = config.pdb_dir(pdb_id)
    query_fasta = path / f"{pdb_id}.fasta"
    for candidate_fasta in config.mer_dir(pdb_id).glob("*.fasta"):
        build_model(pdb_id, query_fasta, candidate_fasta)


def build_model(pdb_id: str, query_fasta: Path, candidate_fasta: Path) -> Path:
    """
    Align the query sequence to one candidate and submit the alignment to
    SWISS-MODEL, downloading the resulting model as
    ``<pdb_id>_<candidate>model.pdb``.
    """
    path = config.pdb_dir(pdb_id)
    candidate_id = candidate_fasta.stem
    aln_path = alignment.write_pair_fasta(query_fasta, candidate_fasta, path / "model.fasta")
    aln_path = alignment.run_clustalw(aln_path)

    output_path = path / f"{pdb_id}_{candidate_id}model.pdb"

    with browser.chrome_driver(headless=False) as bro:
        bro.get(SWISSMODEL_URL)
        bro.find_element(By.ID, "id_sequence_file_upload").send_keys(str(aln_path))
        WebDriverWait(bro, 200).until(
            EC.visibility_of_element_located(
                (By.XPATH, "/html/body/div[2]/div[1]/div[2]/form/div[3]/div/div/div[2]/div/button[2]")
            )
        ).click()

        sleep(INITIAL_WAIT_SECONDS)
        while True:
            try:
                bro.find_element(
                    By.XPATH, '//*[@id="mdl_left_col_01"]/div[1]/div[1]/div/div[1]/button'
                )
            except Exception:
                sleep(POLL_INTERVAL_SECONDS)
                continue
            bro.get(bro.current_url + "01.pdb?display=1")
            model_text = bro.find_element(By.XPATH, "/html/body/pre").text
            output_path.write_text(model_text)
            break

    return output_path
