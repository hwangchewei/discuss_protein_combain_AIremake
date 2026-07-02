"""
Structure surface analysis.

This module is a refactor of the original top-level ``__init__.py``, whose
single ``biopython`` class mixed live logic with several large blocks of
commented-out, superseded scoring formulas. Behaviourally-live code has been
kept (and renamed for clarity); dead/commented branches have been removed.

The class identifies solvent-exposed residues on a protein chain, groups
them into spatially contiguous surface patches, and scores those patches to
suggest which surface region is most likely to mediate homo-oligomer
formation.
"""

from __future__ import annotations

import logging
from typing import List

import numpy
from Bio.PDB import Selection
from Bio.PDB.NeighborSearch import NeighborSearch
from Bio.PDB.PDBParser import PDBParser
from Bio.PDB.ResidueDepth import get_surface, min_dist
from Bio.PDB.SASA import ShrakeRupley
from Bio.PDB.vectors import calc_dihedral

from .amino_acids import CHARGED, HYDROPHILIC, HYDROPHOBIC, SURFACE_WEIGHT_FUNCTION_2

logger = logging.getLogger(__name__)


def _residue_name(residue) -> str:
    """Three-letter residue name, e.g. ``'ALA'``, from a Bio.PDB Residue."""
    return residue.get_resname()


def _residue_seq_id(residue) -> int:
    """Sequence number of a Bio.PDB Residue (its ``id[1]``)."""
    return residue.get_id()[1]


class SurfaceAnalyzer:
    """
    Solvent-accessible surface analysis for a single chain of a PDB structure.

    :param pdb_path: path to a PDB coordinate file.
    :param chain_name: chain identifier to analyse (default ``'A'``).
    """

    def __init__(self, pdb_path: str, chain_name: str = "A") -> None:
        self.chain_name = chain_name
        self.surface_residue_atom: dict = {}
        self.surface_atom: dict = {}
        self.surface_residue_dict: dict = {}

        parser = PDBParser(QUIET=True)
        self.structure = parser.get_structure("structure", pdb_path)
        self.model = self.structure[0]

        try:
            residues = self.model[chain_name].get_list()
            self.has_chain = True
        except KeyError:
            self.has_chain = False
            return

        # Only keep standard (non-heteroatom) residues.
        self.chain = [r for r in residues if r.get_id()[0] == " "]

        try:
            self.surface = get_surface(self.chain)
            self.has_surface = True
        except Exception:  # MSMS can fail on malformed structures
            logger.warning("Could not compute molecular surface for %s chain %s", pdb_path, chain_name)
            self.has_surface = False
            return

        atoms = Selection.unfold_entities(self.chain, "A")
        self.neighbor_search = NeighborSearch(atoms)
        sasa_calculator = ShrakeRupley()
        sasa_calculator.compute(self.model[chain_name], level="C")

    # -- Basic geometry -----------------------------------------------------

    def residue_depth(self, residue):
        """Return ``(min_depth, closest_atom)`` for a residue's atoms vs. the surface."""
        best_depth = float("inf")
        best_atom = None
        for atom in residue.get_unpacked_list():
            depth = min_dist(atom.get_coord(), self.surface)
            if depth < best_depth:
                best_depth = depth
                best_atom = atom
        return best_depth, best_atom

    def surface_residue(self, distance: float) -> List:
        """Return residues whose closest atom lies within ``distance`` of the surface."""
        surface_atoms, surface_residues = [], []
        for residue in self.chain:
            depth, atom = self.residue_depth(residue)
            if depth <= distance:
                surface_atoms.append(atom)
                surface_residues.append(residue)

        self.surface_residue_atom = dict(zip(surface_residues, surface_atoms))
        self.surface_atom = dict(enumerate(surface_atoms))
        self.surface_residue_dict = dict(enumerate(surface_residues))
        return surface_residues

    def total_r_asa(self, residues) -> float:
        """Sum per-atom SASA values across a collection of residues."""
        total = 0.0
        for residue in residues:
            for atom in residue:
                total += round(atom.sasa, 2)
        return total

    # -- Surface patch discovery ---------------------------------------------

    def _search_radius(self) -> float:
        """Radius (Angstrom) used to group nearby surface residues into patches."""
        total_asa = round(self.total_r_asa(self.surface_residue_dict.values()), 0)
        return ((total_asa / 6) / 3.14) ** 0.5 + 2

    def find_surface(self, distance: float) -> List[List]:
        """
        Group solvent-exposed residues into spatially contiguous surface patches.

        Populates and returns ``self.surface_area_list``, a list of residue
        groups (each itself a list of Bio.PDB Residue objects).
        """
        self.surface_residue(distance)
        radius = self._search_radius()

        self.surface_area_list: List[List] = []
        self.surface_area_id_list: List[List[int]] = []

        for index, atom in self.surface_atom.items():
            neighbours = [
                r for r in self.neighbor_search.search(atom.coord, radius, "R")
                if r in self.surface_residue_dict.values()
            ]
            patch = self.surface_area(neighbours, self.surface_residue_dict[index])
            patch_ids = sorted(_residue_seq_id(r) for r in patch)
            if patch_ids not in self.surface_area_id_list:
                self.surface_area_list.append(patch)
                self.surface_area_id_list.append(patch_ids)
        return self.surface_area_list

    def surface_area(self, neighbours: List, residue) -> List:
        """
        Grow a small, roughly-linear patch of surface residues starting from
        ``residue`` using its two nearest surface neighbours as a local frame,
        then extending along that direction using dihedral-angle checks.
        """
        near = []
        for radius in range(1, 10):
            for candidate in self.neighbor_search.search(
                self.surface_residue_atom[residue].coord, radius, "R"
            ):
                if candidate in neighbours and candidate != residue and candidate not in near:
                    near.append(candidate)
                if len(near) == 2:
                    break
            if len(near) == 2:
                break
        else:
            return [residue]

        p1 = self.surface_residue_atom[residue].get_vector()
        p2 = self.surface_residue_atom[near[0]].get_vector()
        p3 = self.surface_residue_atom[near[1]].get_vector()
        patch = [residue, near[0], near[1]]

        positive_hits, negative_hits = 0, 0
        for candidate in neighbours:
            if candidate in patch or candidate not in self.surface_residue_dict.values():
                continue
            p4 = self.surface_residue_atom[candidate].get_vector()
            angle = numpy.degrees(calc_dihedral(p1, p2, p3, p4))
            if abs(angle) > 165 or abs(angle) < 15:
                patch.append(candidate)
            elif angle > 0:
                positive_hits += 1
            else:
                negative_hits += 1

        if positive_hits <= 5 or negative_hits <= 5:
            return patch
        return [residue, near[0], near[1]]

    # -- Patch scoring --------------------------------------------------------

    def surface_area_amount(self, distance: float = 2, weight_threshold: float = 200):
        """
        Score each surface patch found by :meth:`find_surface` using an
        empirical per-residue weighting (:data:`SURFACE_WEIGHT_FUNCTION_2`)
        and return the highest-scoring patch above ``weight_threshold``.

        :returns: ``(description, best_score)`` where ``description`` names
            the residues of the best-scoring patch, or ``(None, 0)`` if no
            patch cleared the threshold.
        """
        self.find_surface(distance)

        best_score = 0.0
        best_patch_description: str | None = None

        for patch in self.surface_area_list:
            weighted_asa = sum(
                self.total_r_asa([residue]) * SURFACE_WEIGHT_FUNCTION_2[_residue_name(residue)]
                for residue in patch
            )
            weighted_asa = round(weighted_asa, 2)

            if weighted_asa >= weight_threshold:
                labelled_ids = sorted(
                    (_residue_seq_id(r), _residue_name(r)) for r in patch
                )
                description = ",".join(f"{seq_id}{name[:1]}" for seq_id, name in labelled_ids)
                if weighted_asa > best_score:
                    best_score = weighted_asa
                    best_patch_description = description

        if best_patch_description is None:
            logger.info("No surface patch cleared the oligomerisation-likelihood threshold")
        return best_patch_description, best_score

    # -- Hydrophobic patch discovery ------------------------------------------

    def find_hydrophobic_patches(self, residues) -> List[List]:
        """
        For each hydrophobic surface residue, grow an expanding-radius search
        until a hydrophilic residue is encountered, keeping only the roughly
        linear portion of the resulting patch (dihedral angle near 0/180deg).
        """
        radius_range = int(((round(self.total_r_asa(self.surface_residue_dict.values()), 0) / 6) / 3.14) ** 0.5 + 1)
        patches: List[List] = []

        for residue in residues:
            if _residue_name(residue) not in HYDROPHOBIC:
                continue
            atom = self.surface_residue_atom[residue]
            for radius in range(1, radius_range):
                patch, hit_hydrophilic = [], False
                for candidate in self.neighbor_search.search(atom.coord, radius, "R"):
                    if candidate in self.surface_residue_dict.values():
                        patch.append(candidate)
                        if _residue_name(candidate) in HYDROPHILIC:
                            hit_hydrophilic = True
                            break
                if not hit_hydrophilic:
                    continue

                to_remove = []
                if len(patch) > 3:
                    p1, p2, p3 = (self.surface_residue_atom[r].get_vector() for r in patch[:3])
                    for extra in patch[3:]:
                        p4 = self.surface_residue_atom[extra].get_vector()
                        angle = abs(numpy.degrees(calc_dihedral(p1, p2, p3, p4)))
                        if not (angle > 170 or angle < 10):
                            to_remove.append(extra)
                for extra in to_remove:
                    patch.remove(extra)
                patches.append(patch)
                break
        return patches

    def hydrophobic_patch_ids(self, distance: float) -> List[List[int]]:
        """
        Convenience wrapper: find surface residues at ``distance``, group
        them into hydrophobic patches, merge collinear patches, and return
        each patch as a sorted list of residue sequence numbers.
        """
        if not getattr(self, "has_chain", False) or not getattr(self, "has_surface", False):
            return []

        residues = self.surface_residue(distance)
        patches = self.find_hydrophobic_patches(residues)
        merged = self._merge_collinear_patches(patches)

        seen: List[List[int]] = []
        for patch in merged:
            ids = sorted({_residue_seq_id(r) for r in patch})
            if ids not in seen:
                seen.append(ids)
        return seen

    def _merge_collinear_patches(self, patches: List[List]) -> List[List]:
        """Merge patches that are roughly collinear extensions of one another."""
        merged = []
        for patch in patches:
            grown = list(patch)
            if len(grown) > 3:
                p1, p2, p3 = (self.surface_residue_atom[r].get_vector() for r in grown[:3])
                for other in patches:
                    if len(other) < 3:
                        if other:
                            p4 = self.surface_residue_atom[other[0]].get_vector()
                            angle = abs(numpy.degrees(calc_dihedral(p1, p2, p3, p4)))
                            if angle > 170 or angle < 10:
                                grown += other
                    else:
                        p4 = self.surface_residue_atom[other[0]].get_vector()
                        p5 = self.surface_residue_atom[other[len(other) // 2]].get_vector()
                        p6 = self.surface_residue_atom[other[-1]].get_vector()
                        avg_angle = (
                            abs(numpy.degrees(calc_dihedral(p1, p2, p3, p4)))
                            + abs(numpy.degrees(calc_dihedral(p1, p2, p3, p5)))
                            + abs(numpy.degrees(calc_dihedral(p1, p2, p3, p6)))
                        ) / 3
                        if avg_angle > 170 or avg_angle < 10:
                            grown += other
            merged.append(grown)
        return merged
