"""
Aggregate reporting across many previously-analysed interfaces.

These functions parse the accumulated text report written by
``interface_composition.interface_residue_composition`` (one growing file,
blocks separated by a line of ``=`` characters, sub-blocks separated by a
line of ``*`` characters) and tally directional trends: for each interface,
did the candidate ("mer") side gain or lose hydrophobic / hydrophilic /
positive / negative residues relative to the query ("mono") side, both per
interface and, in aggregate, per protein.

Equivalent to the original ``interface_change`` and ``interface_amino_change``.
``protein_face`` (surface analysis over whole chains, independent of the
alignment-based interface detection) also lives here since it produces the
same kind of aggregate text report.
"""

from __future__ import annotations

import csv
import logging
from collections import Counter
from pathlib import Path
from typing import List

from . import config
from .amino_acids import HYDROPHILIC, HYDROPHOBIC, NEGATIVE_CHARGE, POSITIVE_CHARGE, letter_counts

logger = logging.getLogger(__name__)

BLOCK_SEPARATOR = "=" * 61
SUBBLOCK_SEPARATOR = "*" * 40


def _parse_summary_line(line: str) -> dict[str, int]:
    """Parse a `LABEL  bic:N(f)  lic:N(f)  p:N(f)  n:N(f)` summary line into counts."""
    parts = line.split("  ")
    counts = {}
    for part in parts[1:]:
        key, value = part.split(":")
        counts[key.strip()] = int(value.split("(")[0])
    return counts


class DirectionTally:
    """Tracks +/-/= counts across the four composition categories (bic/lic/p/n)."""

    CATEGORIES = ("bic", "lic", "p", "n")

    def __init__(self) -> None:
        self.gained = dict.fromkeys(self.CATEGORIES, 0)
        self.lost = dict.fromkeys(self.CATEGORIES, 0)
        self.unchanged = dict.fromkeys(self.CATEGORIES, 0)
        self.direction_patterns: Counter = Counter()

    def record(self, before: dict[str, int], after: dict[str, int]) -> str:
        """Record one before/after comparison; returns the resulting sign pattern, e.g. '+-0+'."""
        pattern = []
        for category in self.CATEGORIES:
            b, a = before.get(category, 0), after.get(category, 0)
            if a > b:
                self.gained[category] += 1
                pattern.append("+")
            elif a < b:
                self.lost[category] += 1
                pattern.append("-")
            else:
                self.unchanged[category] += 1
                pattern.append("0")
        pattern_str = "".join(pattern)
        self.direction_patterns[pattern_str] += 1
        return pattern_str

    def write_csv_rows(self, writer) -> None:
        writer.writerow([self.gained[c] for c in self.CATEGORIES])
        writer.writerow([-self.lost[c] for c in self.CATEGORIES])
        writer.writerow([self.unchanged[c] for c in self.CATEGORIES])
        for pattern, count in sorted(self.direction_patterns.items(), key=lambda kv: kv[1], reverse=True):
            writer.writerow([pattern, count])


def interface_change_report(
    input_report_path: Path | None = None,
    output_csv_path: Path | None = None,
) -> List[str]:
    """
    Tally hydrophobic/hydrophilic/charge gain-or-loss trends across every
    interface recorded in the accumulated interface report, both per
    interface and per protein (i.e. majority direction across all of a
    protein's interfaces).

    Equivalent to the original ``interface_change``.

    :returns: the list of PDB ids that had at least one recorded interface
        (used by ``main`` in the original script to seed further steps).
    """
    input_report_path = input_report_path or config.PROJECT_ROOT / "interface_type_dimer1.txt"
    output_csv_path = output_csv_path or config.PROJECT_ROOT / "change_type1.csv"

    blocks = input_report_path.read_text().split(BLOCK_SEPARATOR)

    per_interface = DirectionTally()
    per_protein = DirectionTally()
    pdb_ids: List[str] = []
    interface_count = 0
    protein_count = 0

    for block in blocks:
        lines = block.split("\n")
        if len(lines) <= 5:
            continue
        pdb_ids.append(lines[2])
        protein_count += 1
        protein_tally = DirectionTally()

        for sub_block in block.split(SUBBLOCK_SEPARATOR):
            lines = sub_block.split("\n")
            if len(lines) < 7:
                continue
            try:
                before = _parse_summary_line(lines[-4])
                after = _parse_summary_line(lines[-3])
            except (IndexError, ValueError):
                continue
            interface_count += 1
            pattern = protein_tally.record(before, after)
            per_interface.direction_patterns[pattern] += 1

        for category in DirectionTally.CATEGORIES:
            per_interface.gained[category] += protein_tally.gained[category]
            per_interface.lost[category] += protein_tally.lost[category]
            per_interface.unchanged[category] += protein_tally.unchanged[category]

        # A protein's overall direction per category is whichever of
        # gained/lost occurred across more of *its own* interfaces.
        per_protein.record(protein_tally.lost, protein_tally.gained)

    with output_csv_path.open("w", newline="") as fp:
        writer = csv.writer(fp)
        _write_full(writer, per_interface, interface_count, per_protein, protein_count)

    return pdb_ids


def _write_full(writer, per_interface: DirectionTally, interface_count: int, per_protein: DirectionTally, protein_count: int) -> None:
    per_interface.write_csv_rows(writer)
    writer.writerow([interface_count])
    per_protein.write_csv_rows(writer)
    writer.writerow([protein_count])


def interface_amino_change_report(
    input_report_path: Path | None = None,
    output_csv_path: Path | None = None,
) -> None:
    """
    Tally which specific amino-acid substitutions occur between the query
    and candidate side of each interface residue pair, plus overall
    per-letter frequency and composition-class counts.

    Equivalent to the original ``interface_amino_change``.
    """
    input_report_path = input_report_path or config.PROJECT_ROOT / "interface_type_dimer1.txt"
    output_csv_path = output_csv_path or config.PROJECT_ROOT / "amino_change_type1.csv"

    blocks = input_report_path.read_text().split(BLOCK_SEPARATOR)

    substitution_counts: Counter = Counter()
    substitution_type_counts: Counter = Counter()
    same_pairs: List[List[str]] = [["same pdbid"]]
    query_letters: Counter = Counter()
    candidate_letters: Counter = Counter()
    query_classes = dict(bic=0, lic=0, p=0, n=0, gap=0, other=0)
    candidate_classes = dict(bic=0, lic=0, p=0, n=0, gap=0, other=0)

    for block in blocks:
        if len(block.split("\n")) <= 5:
            continue
        for sub_block in block.split(SUBBLOCK_SEPARATOR):
            lines = sub_block.split("\n")
            if len(lines) < 7:
                continue
            query_residues = _parse_residue_list(lines[-7])
            candidate_residues = _parse_residue_list(lines[-6])

            if query_residues == candidate_residues:
                same_pairs.append([lines[-4].split("  ")[0], lines[-3].split("  ")[0]])
                substitution_type_counts["same"] += 1
                continue

            for query_res, candidate_res in zip(query_residues, candidate_residues):
                if query_res == candidate_res:
                    continue
                substitution_counts[f"{query_res} -> {candidate_res}"] += 1
                _tally_letter(query_res, query_letters, query_classes)
                _tally_letter(candidate_res, candidate_letters, candidate_classes)
                substitution_type_counts[
                    f"{_class_label(query_res)} -> {_class_label(candidate_res)}"
                ] += 1

    with output_csv_path.open("w", newline="") as fp:
        writer = csv.writer(fp)
        for letter in "ACDEFGHKILMNPQRSTVYW":
            writer.writerow([letter, query_letters.get(letter, 0), candidate_letters.get(letter, 0)])
        for substitution, count in sorted(substitution_counts.items(), key=lambda kv: kv[1], reverse=True):
            writer.writerow([substitution, count])
        for sub_type, count in sorted(substitution_type_counts.items(), key=lambda kv: kv[1], reverse=True):
            writer.writerow([sub_type, count])
        writer.writerows(same_pairs)
        writer.writerow([query_classes[k] for k in ("bic", "lic", "p", "n", "gap", "other")])
        writer.writerow([candidate_classes[k] for k in ("bic", "lic", "p", "n", "gap", "other")])


def _parse_residue_list(line: str) -> List[str]:
    return line.replace("[", "").replace("]", "").replace("'", "").replace(" ", "").split(",")


def _class_label(one_letter_residue: str) -> str:
    from .amino_acids import ONE_TO_THREE, classify_residue

    if one_letter_residue == "-":
        return "-"
    three = ONE_TO_THREE.get(one_letter_residue)
    if three is None:
        return one_letter_residue
    return classify_residue(three).value


def protein_face_report(pdb_id: str, output_path: Path | None = None, surface_distance: float = 4.5) -> None:
    """
    Report the hydrophobic/hydrophilic/positive/negative composition of the
    solvent-exposed surface of every single-chain PDB file in ``mer/``.

    Equivalent to the original ``protein_face``. Uses
    :class:`structure_surface.SurfaceAnalyzer` rather than the ``biopython``
    class from the original top-level ``__init__.py``.
    """
    if len(pdb_id) != 4:
        return

    from .amino_acids import composition
    from .structure_surface import SurfaceAnalyzer

    output_path = output_path or config.PROJECT_ROOT / "protein_face_type.txt"
    mer_directory = config.mer_dir(pdb_id)

    ids: List[str] = []
    summary_lines: List[str] = []

    for pdb_path in mer_directory.glob("*.pdb"):
        if "_" in pdb_path.stem:
            continue
        analyzer = SurfaceAnalyzer(str(pdb_path))
        if not getattr(analyzer, "has_chain", False) or not getattr(analyzer, "has_surface", False):
            continue
        surface_residues = analyzer.surface_residue(surface_distance)
        counts = composition([r.get_resname() for r in surface_residues])
        total = sum(counts.values()) or 1

        def part(class_value: str) -> str:
            n = next((v for k, v in counts.items() if k.value == class_value), 0)
            return f"{round(n / total, 3)}({n})"

        ids.append(pdb_path.stem + " ")
        summary_lines.append(
            "{:5s}{:^16s}{:^16s}{:^16s}{:^16s}".format(
                pdb_path.stem, part("Bic"), part("Lic"), part("P"), part("N")
            )
        )

    with output_path.open("a") as fp:
        fp.write(pdb_id + "\n")
        fp.writelines(ids)
        fp.write("\n")
        fp.write("{:5s}{:^16s}{:^16s}{:^16s}{:^16s}\n".format("ID", "hydrophobic", "hydrophilic", "positive", "negative"))
        fp.write("\n".join(summary_lines) + "\n")
        fp.write("*" * 70 + "\n\n")


def _tally_letter(one_letter_residue: str, letter_counter: Counter, class_counts: dict[str, int]) -> None:
    if one_letter_residue == "-":
        class_counts["gap"] += 1
        return
    letter_counter[one_letter_residue] += 1
    from .amino_acids import ONE_TO_THREE

    three = ONE_TO_THREE.get(one_letter_residue)
    if three in HYDROPHOBIC:
        class_counts["bic"] += 1
    elif three in HYDROPHILIC:
        class_counts["lic"] += 1
    elif three in POSITIVE_CHARGE:
        class_counts["p"] += 1
    elif three in NEGATIVE_CHARGE:
        class_counts["n"] += 1
    else:
        class_counts["other"] += 1
