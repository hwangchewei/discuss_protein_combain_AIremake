"""
EBI PISA interface analysis.

PISA (Proteins, Interfaces, Structures and Assemblies) reports, for a PDB
entry, the geometric/energetic properties of every pairwise molecule-molecule
interface in the asymmetric unit. This module:

- downloads PISA's pre-computed interface XML for an entry
  (:func:`download_interfaces_xml`, was ``pisa_download``);
- filters out interfaces that aren't genuine protein-protein contacts, or
  that involve a ligand (:func:`filter_true_protein_interfaces`, was
  ``check_mer``);
- strips ligand/heteroatom records from PDB coordinate files
  (:func:`strip_ligands`, was ``pdb_del_ligand``);
- drives the interactive PISA web tool for pairs it doesn't have
  pre-computed data for (:func:`analyse_pair_interactively`, was ``pisa``),
  and checks whether a pair is *already* a known assembly
  (:func:`check_already_known_assembly`, was ``pisa_mer_check`` /
  ``pisa_mer_check_``).

The interactive-scraping functions are inherently coupled to PISA's current
HTML structure (hence the long XPath strings); they're kept close to the
original logic but with duplicated file-moving code factored out.
"""

from __future__ import annotations

import csv
import logging
import shutil
from pathlib import Path
from time import sleep

from lxml import etree
from selenium.webdriver.common.by import By
from selenium.webdriver.support import expected_conditions as EC
from selenium.webdriver.support.ui import WebDriverWait

from . import browser, config, http_client

logger = logging.getLogger(__name__)

MIN_VALID_PDB_SIZE_BYTES = 10_000


def download_interfaces_xml(pdb_id: str) -> None:
    """Download PISA's interface-analysis XML for every candidate in ``<pdb_id>/``."""
    if len(pdb_id) != 4:
        return
    path = config.pdb_dir(pdb_id)
    for entry in path.iterdir():
        if entry.name in ("no_aligne", "monomer"):
            continue
        candidate_id = entry.name.split("_")[0]
        http_client.download_to_file(
            http_client.fetch_pisa_interfaces_xml, candidate_id, path / f"{candidate_id}.xml"
        )


def _move_related_files(mer_dir: Path, stem: str, destination: Path) -> None:
    """Move every ``<stem>.*`` and ``<parent>_<stem>.pdb`` file into ``destination``."""
    destination.mkdir(parents=True, exist_ok=True)
    for pattern in (f"{stem}.xml", f"{stem}.fasta", f"{stem}.pdb", f"*_{stem}.pdb"):
        for match in mer_dir.glob(pattern):
            shutil.move(str(match), str(destination))


def filter_true_protein_interfaces(pdb_id: str) -> None:
    """
    Discard candidates whose PISA XML shows no genuine, non-ligand
    protein-protein interface (crystallographic-symmetry interfaces below
    the solvation-energy threshold, or interfaces mediated by a ligand).

    Equivalent to the original ``check_mer``.
    """
    if pdb_id in ("oligState", "no_homomer", "no", "no_aligne", "monomer", "false_mer") or "." in pdb_id:
        return

    mer_directory = config.mer_dir(pdb_id)
    false_mer_dir = mer_directory / "false_mer"
    false_mer_dir.mkdir(parents=True, exist_ok=True)

    any_xml = False
    for xml_path in mer_directory.glob("*.xml"):
        any_xml = True
        stem = xml_path.stem
        tree = etree.parse(str(xml_path))
        interfaces = tree.xpath("//interface")

        is_true_interface = False
        is_ligand_mediated = False
        ligand_chains = []

        for i in range(1, len(interfaces) + 1):
            solvation_energy = float(tree.xpath(f"//interface[{i}]//int_solv_en/text()")[0])
            symop1 = tree.xpath(f"//interface[{i}]//molecule[1]/symop/text()")[0]
            symop2 = tree.xpath(f"//interface[{i}]//molecule[2]/symop/text()")[0]
            class1 = tree.xpath(f"//interface[{i}]//molecule[1]/class/text()")[0]
            class2 = tree.xpath(f"//interface[{i}]//molecule[2]/class/text()")[0]
            chain1 = tree.xpath(f"//interface[{i}]//molecule[1]/chain_id/text()")[0]
            chain2 = tree.xpath(f"//interface[{i}]//molecule[2]/chain_id/text()")[0]

            same_asymmetric_unit = symop1 == "x,y,z" and symop2.lower() == "x,y,z"

            if same_asymmetric_unit and class1 == class2 == "Protein" and solvation_energy < config.MAX_INTERFACE_SOLVATION_ENERGY:
                is_true_interface = True
                break

            if same_asymmetric_unit and ("Ligand" in (class1, class2)) and class1 != class2:
                is_ligand_mediated = True
                ligand_chains.append(chain1 if class1 == "Ligand" else chain2)

        if ligand_chains:
            (mer_directory / f"{stem}_ligand.txt").write_text(
                stem + "\n" + " ".join(str(c) for c in ligand_chains) + "\n"
            )

        if is_ligand_mediated:
            logger.info("%s/%s: ligand-mediated interface, moving to false_mer", pdb_id, stem)
            _move_related_files(mer_directory, stem, false_mer_dir)

    if not any_xml:
        shutil.move(str(config.pdb_dir(pdb_id)), str(config.bucket_dir("false_mer")))


def strip_ligands(pdb_id: str) -> None:
    """
    Write a ``<name>_del.pdb`` copy of every plain (non-underscore-suffixed)
    PDB file in ``<pdb_id>/mer/`` keeping only ATOM/CRYST1/TER/END records.

    Equivalent to the original ``pdb_del_ligand``.
    """
    if len(pdb_id) != 4:
        return
    mer_directory = config.mer_dir(pdb_id)
    keep_prefixes = ("ATOM", "CRYST1", "TER", "END")
    for pdb_path in mer_directory.glob("*.pdb"):
        if "_" in pdb_path.stem:
            continue
        lines = pdb_path.read_text().split("\n")
        kept = [line for line in lines if any(line.startswith(p) or p in line for p in keep_prefixes)]
        (mer_directory / f"{pdb_path.stem}_del.pdb").write_text("\n".join(kept) + "\n")


def move_promoted_files(pdb_id: str) -> None:
    """
    Promote files that were sorted into ``mer/alreadymer/`` back up into
    ``mer/``, then clean up the now-empty ``alreadymer`` directory.

    Equivalent to the original ``move_file``.
    """
    if len(pdb_id) != 4:
        return
    mer_directory = config.mer_dir(pdb_id)
    already_mer_dir = mer_directory / "alreadymer"
    if not already_mer_dir.is_dir():
        return
    for entry in list(already_mer_dir.iterdir()):
        shutil.move(str(entry), str(mer_directory))
    if not any(already_mer_dir.iterdir()):
        already_mer_dir.rmdir()


# ---------------------------------------------------------------------------
# Interactive PISA-web checks (require a browser)
# ---------------------------------------------------------------------------

def check_already_known_assembly(pdb_id: str) -> None:
    """
    For every ``<pdb_id>_<candidate>.pdb`` (not yet ``_del``) in ``mer/``,
    upload it to the PISA web tool to check whether it already matches a
    known crystallographic assembly, then bucket the whole entry into
    ``alreadymer`` if nothing is left needing further analysis.

    Equivalent to the original ``pisa_mer_check`` (dispatch loop).
    """
    if len(pdb_id) != 4:
        return
    mer_directory = config.mer_dir(pdb_id)
    for pdb_path in mer_directory.glob("*.pdb"):
        if "_" in pdb_path.stem and "del" not in pdb_path.stem:
            _check_pair_is_already_assembly(pdb_id, pdb_path.name)

    if not any(mer_directory.glob("*.pdb")):
        shutil.move(str(config.pdb_dir(pdb_id)), str(config.bucket_dir("alreadymer")))


def _check_pair_is_already_assembly(pdb_id: str, filename: str, error_count: int = 0) -> None:
    """Equivalent to the original ``pisa_mer_check_``."""
    mer_directory = config.mer_dir(pdb_id)
    already_mer_dir = mer_directory / "alreadymer"
    already_mer_dir.mkdir(parents=True, exist_ok=True)

    with browser.chrome_driver(headless=True) as bro:
        bro.get("https://www.ebi.ac.uk/pdbe/pisa/")
        try:
            WebDriverWait(bro, 200).until(
                EC.visibility_of_element_located(
                    (By.XPATH, "/html/body/div[2]/div[2]/div/table/tbody/tr/td[2]/table/tbody/tr/td/div/div/div/form/span/span/button")
                )
            ).click()
            WebDriverWait(bro, 1).until(
                EC.visibility_of_element_located((By.XPATH, '//*[@id="sform"]/tbody/tr[4]/td/u/input'))
            ).click()
            WebDriverWait(bro, 10).until(
                EC.visibility_of_element_located((By.XPATH, '//*[@id="sform"]/tbody/tr[3]/td/b/input[2]'))
            ).send_keys(str(mer_directory / filename))
            WebDriverWait(bro, 10).until(
                EC.visibility_of_element_located((By.XPATH, '//*[@id="sform"]/tbody/tr[3]/td/b/input[3]'))
            ).click()
            WebDriverWait(bro, 300).until(
                EC.visibility_of_element_located((By.XPATH, "/html/body/div[2]/div[2]/div/form/table/tbody//td/span[3]/span/button"))
            ).click()
            WebDriverWait(bro, 200).until(
                EC.visibility_of_element_located((By.XPATH, "/html/body/div[2]/div[2]/div/form/table[2]/tbody"))
            )
        except Exception:
            if error_count > 3:
                config.PROJECT_ROOT.joinpath("alreadymer_timeout.txt").open("a").write(f"{pdb_id} {filename}\n")
                return
            return _check_pair_is_already_assembly(pdb_id, filename, error_count + 1)

        rows = bro.find_elements(By.XPATH, "/html/body/div[2]/div[2]/div/form/table[2]/tbody/tr")
        found_assembly = False
        csv_path = mer_directory / f"{pdb_id}_alreadymer.csv"

        for row_index in range(len(rows)):
            sub_rows = bro.find_elements(
                By.XPATH, f"/html/body/div[2]/div[2]/div/form/table[2]/tbody/tr[{row_index}]//tr"
            )
            if not sub_rows:
                continue
            found_assembly = True
            with csv_path.open("a", newline="") as fp:
                csv.writer(fp).writerow([filename.split("_")[1].split(".")[0]])
            for sub_row in range(3, len(sub_rows) + 1):
                text = bro.find_element(
                    By.XPATH,
                    f"/html/body/div[2]/div[2]/div/form/table[2]/tbody/tr[{row_index}]//tr[{sub_row}]",
                ).text
                with csv_path.open("a", newline="") as fp:
                    row = text.split()
                    csv.writer(fp).writerow(([" "] + row) if sub_row > 3 else row)

        if found_assembly:
            candidate_id = filename.split(".")[0].split("_")[1]
            for match in mer_directory.glob(f"*{candidate_id}*"):
                shutil.move(str(match), str(already_mer_dir))


def analyse_pair_interactively(pdb_id: str) -> None:
    """
    Drive the interactive PISA web tool for every ``<pdb_id>_<candidate>.pdb``
    that doesn't already have a cached ``.csv`` result, recording per-interface
    solvation energy / hydrogen-bond / salt-bridge / disulphide counts.

    Equivalent to the original ``pisa`` function. This is the most UI-coupled
    part of the pipeline (PISA exposes no stable API for this step); expect
    to need to update the XPath expressions if PISA's page layout changes.
    """
    if len(pdb_id) != 4:
        return

    mer_directory = config.mer_dir(pdb_id)
    if any(mer_directory.glob("*.csv")):
        return

    for stale_xml in mer_directory.glob("*_*.xml"):
        stale_xml.unlink()

    csv_path = mer_directory / f"{pdb_id}.csv"
    if csv_path.exists():
        csv_path.unlink()

    false_mer_dir = mer_directory / "false_mer"
    monomer_dir = mer_directory / "monomer"
    overlapping_dir = mer_directory / "overlapping"

    for pdb_path in list(mer_directory.glob("*_*.pdb")):
        candidate_id = pdb_path.stem.split("_")[1]

        if pdb_path.stat().st_size < MIN_VALID_PDB_SIZE_BYTES:
            _move_related_files(mer_directory, candidate_id, false_mer_dir)
            shutil.move(str(pdb_path), str(false_mer_dir))
            continue

        with browser.chrome_driver(headless=True) as bro:
            bro.get("https://www.ebi.ac.uk/pdbe/pisa/")
            WebDriverWait(bro, 200).until(
                EC.visibility_of_element_located(
                    (By.XPATH, "/html/body/div[2]/div[2]/div/table/tbody/tr/td[2]/table/tbody/tr/td/div/div/div/form/span/span/button")
                )
            ).click()
            WebDriverWait(bro, 1).until(
                EC.visibility_of_element_located((By.XPATH, '//*[@id="sform"]/tbody/tr[4]/td/u/input'))
            ).click()
            WebDriverWait(bro, 200).until(
                EC.visibility_of_element_located((By.XPATH, '//*[@id="sform"]/tbody/tr[3]/td/b/input[2]'))
            ).send_keys(str(pdb_path))
            WebDriverWait(bro, 200).until(
                EC.visibility_of_element_located((By.XPATH, '//*[@id="sform"]/tbody/tr[3]/td/b/input[3]'))
            ).click()
            WebDriverWait(bro, 200).until(
                EC.visibility_of_element_located((By.XPATH, '//*[@id="pdbeSubmitButton"]/span/button'))
            ).click()

            try:
                sleep(60)
                bro.find_element(By.XPATH, "/html/body/div[2]/div[2]/div/div/form/input[5]")
                monomer_dir.mkdir(parents=True, exist_ok=True)
                shutil.move(str(pdb_path), str(monomer_dir))
                _move_related_files(mer_directory, candidate_id, monomer_dir)
                continue
            except Exception:
                pass

            WebDriverWait(bro, 600).until(
                EC.visibility_of_element_located((By.XPATH, '//*[@id="makePageHeader0"]//tbody/tr[4]/td[2]/span/span/button'))
            ).click()

            overlapping, overlapping_count = "No overlapping", "0"
            try:
                WebDriverWait(bro, 3).until(
                    EC.visibility_of_element_located(
                        (By.XPATH, "/html/body/div[2]/div[2]/div/form/table[2]/tbody/tr[3]/td/table/tbody/tr[2]/td/font/font/p[1]/strong")
                    )
                )
                overlapping = bro.find_element(
                    By.XPATH, "/html/body/div[2]/div[2]/div/form/table[2]/tbody/tr[3]/td/table/tbody/tr[2]/td/font/font/p[1]/strong"
                ).text
                if overlapping == "Overlapping structures":
                    overlapping_dir.mkdir(parents=True, exist_ok=True)
                    shutil.move(str(pdb_path), str(overlapping_dir))
                    _move_related_files(mer_directory, candidate_id, overlapping_dir)
                    bro.find_element(By.XPATH, '//*[@id="makePageHeader0"]//tbody/tr[2]/td[2]/span/span/button').click()
                    overlapping_count = bro.find_element(
                        By.XPATH, "/html/body/div[2]/div[2]/div/form/table[2]/tbody/tr/td[2]"
                    ).text
                    bro.find_element(By.XPATH, '//*[@id="makePageHeader0"]//tbody/tr[2]/td[3]/span/span/button').click()
            except Exception:
                bro.find_element(By.XPATH, '//*[@id="makePageHeader0"]//tbody/tr[2]/td[2]/span/span/button').click()

            try:
                WebDriverWait(bro, 300).until(
                    EC.visibility_of_element_located((By.XPATH, "/html/body/div[2]/div[2]/div/form/table[2]/tbody/tr[1]/td/span[3]/span/button"))
                )
            except Exception:
                continue  # retry this candidate on a future run

            with csv_path.open("a", newline="") as fp:
                csv.writer(fp).writerow([candidate_id, overlapping, overlapping_count])

            _scrape_interfaces(bro, mer_directory, pdb_id, candidate_id, monomer_dir)

    if not any(mer_directory.glob("*.pdb")):
        shutil.move(str(config.pdb_dir(pdb_id)), str(config.bucket_dir("overlapping")))


def _scrape_interfaces(bro, mer_directory: Path, pdb_id: str, candidate_id: str, monomer_dir: Path) -> None:
    """Scrape each listed interface row of a PISA result page, saving true (non-symmetry) contacts."""
    interface_rows = bro.find_elements(
        By.XPATH, "/html/body/div[2]/div[2]/div/form/table[2]/tbody/tr[2]/td/table/tbody/tr"
    )
    found_true_interface = False

    for row_index in range(1, len(interface_rows) + 1):
        cells = bro.find_elements(
            By.XPATH, f'//*[@id="content"]/div/form/table[2]/tbody/tr[2]/td/table/tbody/tr[{row_index}]/td'
        )
        # PISA's column layout shifts depending on whether a leading
        # checkbox column is present (20 vs 21 total <td> cells).
        if len(cells) == 20:
            col = {"n": 1, "p1": 3, "p2": 8, "symop": 9, "hbonds": 17, "saltbridges": 18, "ssbonds": 19}
        elif len(cells) == 21:
            col = {"n": 2, "p1": 4, "p2": 9, "symop": 10, "hbonds": 18, "saltbridges": 19, "ssbonds": 20}
        else:
            continue

        def cell(name: str) -> str:
            return bro.find_element(
                By.XPATH,
                f'//*[@id="content"]/div/form/table[2]/tbody/tr[2]/td/table/tbody/tr[{row_index}]/td[{col[name]}]',
            ).text

        p1, p2, symop = cell("p1"), cell("p2"), cell("symop")
        if ":" in p1 or ":" in p2 or symop != "  x,y,z  ":
            continue

        found_true_interface = True
        bro.find_element(
            By.XPATH, f'//*[@id="content"]/div/form/table[2]/tbody/tr[2]/td/table/tbody/tr[{row_index}]/td[{col["n"]}]'
        ).click()
        WebDriverWait(bro, 200).until(
            EC.visibility_of_element_located((By.XPATH, "/html/body/div[2]/div[2]/div/form/table[2]/tbody/tr/td/span[1]/span/button"))
        ).click()
        bro.switch_to.window(bro.window_handles[-1])
        interface_xml = bro.find_element(By.CLASS_NAME, "pretty-print").text
        (mer_directory / f"{candidate_id}_{p1}{p2}.xml").write_text(interface_xml)
        with (mer_directory / f"{pdb_id}.csv").open("a", newline="") as fp:
            csv.writer(fp).writerow(
                [cell("n"), p1, p2, symop, cell("hbonds"), cell("saltbridges"), cell("ssbonds")]
            )
        bro.close()
        bro.switch_to.window(bro.window_handles[0])
        bro.back()

    if not found_true_interface:
        monomer_dir.mkdir(parents=True, exist_ok=True)
        _move_related_files(mer_directory, candidate_id, monomer_dir)
