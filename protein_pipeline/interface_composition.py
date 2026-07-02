"""
Amino-acid composition analysis for candidate interfaces.

Two levels of analysis from the original script live here:

- :func:`whole_sequence_composition` (was ``protein_type``): coarse
  hydrophobic/hydrophilic fraction of each *whole* candidate sequence.
- :func:`interface_residue_composition` (was ``interface_check``): a much
  more involved routine that aligns the query sequence to each mer
  candidate, cross-references PISA's per-residue buried-surface-area (BSA)
  data to find which residues actually sit at the interface, and reports
  the hydrophobic / hydrophilic / positive / negative composition of that
  interface on both sides of the pair.

The residue-numbering reconciliation in ``interface_check`` (PISA numbers
residues by their position in the deposited PDB SEQRES/ATOM records, which
can be offset from the ClustalW alignment's numbering by a fixed but a
priori unknown shift) is preserved as :func:`_resolve_numbering_offset`.
"""

from __future__ import annotations

import logging
from dataclasses import dataclass, field
from pathlib import Path
from typing import List

from lxml import etree

from . import alignment, config
from .amino_acids import ONE_TO_THREE, THREE_TO_ONE, classify_residue, composition

logger = logging.getLogger(__name__)


def whole_sequence_composition(pdb_id: str, output_path: Path | None = None) -> None:
    """
    Report the hydrophobic/hydrophilic fraction of the query sequence and of
    each ``mer/*.fasta`` candidate sequence. Equivalent to the original
    ``protein_type``.
    """
    if len(pdb_id) != 4:
        return

    output_path = output_path or config.PROJECT_ROOT / "protein_type_new.txt"
    path = config.pdb_dir(pdb_id)
    query_sequence = (path / f"{pdb_id}.fasta").read_text().split("\n")[1]

    lines = [_format_fraction_line(pdb_id, query_sequence)]
    ids = [pdb_id]

    for fasta_path in config.mer_dir(pdb_id).glob("*.fasta"):
        try:
            sequence = fasta_path.read_text().split("\n")[1]
        except IndexError:
            continue
        ids.append(fasta_path.stem)
        lines.append(_format_fraction_line(fasta_path.stem, sequence))

    with output_path.open("a") as fp:
        fp.write(pdb_id + "\n")
        fp.writelines(f"{i} " for i in ids)
        fp.write("\n")
        fp.write("{:5s}{:^16s}{:^16s}\n".format("ID", "hydrophobic", "hydrophilic"))
        fp.write("\n".join(lines) + "\n")
        fp.write("*" * 70 + "\n\n")


def _format_fraction_line(label: str, sequence: str) -> str:
    hydrophobic = hydrophilic = 0
    for one_letter in sequence:
        three_letter = ONE_TO_THREE.get(one_letter)
        if three_letter is None:
            continue
        residue_class = classify_residue(three_letter)
        if residue_class.value == "Bic":
            hydrophobic += 1
        elif residue_class.value in ("Lic", "P", "N"):
            hydrophilic += 1
    length = len(sequence) or 1
    return "{:5s}{:^16.3f}{:^16.3f}".format(label, round(hydrophobic / length, 3), round(hydrophilic / length, 3))


@dataclass
class InterfaceCompositionResult:
    """Composition summary for one query/candidate interface pair."""

    candidate_id: str
    chain_query: str
    chain_candidate: str
    query_residues: List[str]
    candidate_residues: List[str]
    residue_numbers: List[int]
    unmatched: List[str] = field(default_factory=list)


def _find_best_true_interface(tree) -> int | None:
    """Return the 1-based index of the lowest-solvation-energy true protein interface, if any."""
    interfaces = tree.xpath("//interface")
    scored = {}
    for i in range(1, len(interfaces) + 1):
        solvation_energy = float(tree.xpath(f"//interface[{i}]//int_solv_en/text()")[0])
        symop1 = tree.xpath(f"//interface[{i}]//molecule[1]/symop/text()")[0]
        symop2 = tree.xpath(f"//interface[{i}]//molecule[2]/symop/text()")[0]
        class1 = tree.xpath(f"//interface[{i}]//molecule[1]/class/text()")[0]
        class2 = tree.xpath(f"//interface[{i}]//molecule[2]/class/text()")[0]
        if symop1 == "x,y,z" and symop2.lower() == "x,y,z" and class1 == class2 == "Protein":
            scored[solvation_energy] = i
    if not scored:
        return None
    return scored[min(scored)]


def _resolve_numbering_offset(pisa_residue_names: List[str], mer_index: dict[int, str], start_hint: int = 1, window: int = 100) -> int:
    """
    PISA's residue numbering and the ClustalW-alignment-derived numbering
    (``mer_index``) can differ by a constant offset. Find that offset by
    sliding a 3-residue window of the PISA residue list against
    ``mer_index`` until three consecutive residues match.
    """
    for offset in range(start_hint, start_hint + window):
        try:
            if (
                THREE_TO_ONE.get(pisa_residue_names[0]) == mer_index.get(offset)
                and THREE_TO_ONE.get(pisa_residue_names[1]) == mer_index.get(offset + 1)
                and THREE_TO_ONE.get(pisa_residue_names[2]) == mer_index.get(offset + 2)
            ):
                return offset
        except IndexError:
            break
    return start_hint


def interface_residue_composition(pdb_id: str, report_path: Path | None = None) -> None:
    """
    For each aligned mer candidate under ``mer/dimer/really_dimer/``, align
    it to the query sequence, map PISA's buried-residue list onto that
    alignment, and report the hydrophobic/hydrophilic/charged composition of
    the actual interface residues on both sides.

    Equivalent to the original ``interface_check``. Results are appended to
    ``report_path`` (default: ``<PROJECT_ROOT>/interface_type_test.txt``),
    matching the accumulation-by-append behaviour ``main()`` relied on.
    """
    if len(pdb_id) != 4:
        return

    report_path = report_path or config.PROJECT_ROOT / "interface_type_test.txt"
    candidates_dir = config.pdb_dir(pdb_id) / "mer" / "dimer" / "really_dimer"
    if not candidates_dir.is_dir():
        return
    candidates = list(candidates_dir.glob("*.fasta"))
    if not candidates:
        return

    error_dir = candidates_dir / "interface_error"
    error_dir.mkdir(parents=True, exist_ok=True)
    fasta_workdir = candidates_dir / "fasta"
    fasta_workdir.mkdir(parents=True, exist_ok=True)

    query_fasta = config.pdb_dir(pdb_id) / f"{pdb_id}.fasta"

    for candidate_fasta in candidates:
        candidate_id = candidate_fasta.stem
        try:
            _process_one_interface(pdb_id, candidate_id, query_fasta, candidate_fasta, candidates_dir, fasta_workdir, report_path)
        except Exception:
            logger.exception("Failed processing interface %s / %s", pdb_id, candidate_id)


def _process_one_interface(
    pdb_id: str,
    candidate_id: str,
    query_fasta: Path,
    candidate_fasta: Path,
    candidates_dir: Path,
    fasta_workdir: Path,
    report_path: Path,
) -> None:
    pairwise = alignment.align_pair(query_fasta, candidate_fasta, fasta_workdir, f"{pdb_id}_{candidate_id}")
    query_index, mer_index = pairwise.residue_index_map()

    xml_path = candidates_dir / f"{candidate_id}.xml"
    tree = etree.parse(str(xml_path))
    interface_number = _find_best_true_interface(tree)

    no_contact_dir = candidates_dir / "XYZnocontact"
    if interface_number is None:
        no_contact_dir.mkdir(parents=True, exist_ok=True)
        for related in candidates_dir.glob(f"*{candidate_id}*"):
            related.rename(no_contact_dir / related.name)
        return

    residue_numbers = tree.xpath(f"//interface[{interface_number}]//molecule[1]/residues/residue/seq_num/text()")
    residue_names = tree.xpath(f"//interface[{interface_number}]//molecule[1]/residues/residue/name/text()")
    buried_areas = tree.xpath(f"//interface[{interface_number}]//molecule[1]//residues/residue/bsa/text()")
    chain_query = tree.xpath(f"//interface[{interface_number}]//molecule[1]/chain_id/text()")[0]
    chain_candidate = tree.xpath(f"//interface[{interface_number}]//molecule[2]/chain_id/text()")[0]

    offset = _resolve_numbering_offset(residue_names, mer_index)

    query_interface_residues, candidate_interface_residues, seq_numbers, unmatched = [], [], [], []

    for i, residue_name in enumerate(residue_names):
        position = offset + i
        if position > len(mer_index):
            break
        if int(residue_numbers[i]) == 0:
            continue
        if residue_name not in THREE_TO_ONE:
            if float(buried_areas[i]) != 0:
                unmatched.append(f"{residue_name} {position} {query_index.get(position, '?')}")
            continue

        one_letter = THREE_TO_ONE[residue_name]
        if one_letter != mer_index.get(position):
            offset = _resolve_numbering_offset(residue_names[i:], mer_index, start_hint=offset + i) - i

        if float(buried_areas[i]) == 0:
            continue

        position = offset + i
        if THREE_TO_ONE.get(residue_name) == mer_index.get(position):
            query_interface_residues.append(query_index[position])
            candidate_interface_residues.append(mer_index[position])
            seq_numbers.append(position)
        else:
            unmatched.append(f"{residue_name} {position}")

    query_summary = _format_composition_summary(pdb_id, query_interface_residues, one_letter=True)
    candidate_summary = _format_composition_summary(candidate_id, candidate_interface_residues, one_letter=True)

    with report_path.open("a") as fp:
        fp.write(f"{candidate_id} {chain_query} {chain_candidate}\n")
        fp.write(pairwise.query_aligned + "\n")
        fp.write(pairwise.mer_aligned + "\n")
        fp.write(str(query_interface_residues) + "\n")
        fp.write(str(candidate_interface_residues) + "\n")
        fp.write(str(seq_numbers) + "\n")
        fp.write(query_summary + "\n")
        fp.write(candidate_summary + "\n")
        fp.write(str(unmatched) + "\n")
        fp.write("*" * 40 + "\n")


def _format_composition_summary(label: str, residues: List[str], *, one_letter: bool) -> str:
    """Render a ``LABEL  bic:N(f)  lic:N(f)  p:N(f)  n:N(f)`` summary line."""
    valid = [r for r in residues if r != "-" and r.lower() != "x"]
    counts = composition(valid, one_letter=one_letter)
    total = sum(counts.values()) or 1

    def part(class_value: str) -> str:
        n = next((v for k, v in counts.items() if k.value == class_value), 0)
        return f"{n}({round(n / total, 3)})"

    return "{}  bic:{}  lic:{}  p:{}  n:{}".format(label, part("Bic"), part("Lic"), part("P"), part("N"))
