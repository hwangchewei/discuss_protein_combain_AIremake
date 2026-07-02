"""
Pairwise sequence alignment via ClustalW.

The original script called ``ClustalwCommandline`` and then hand-parsed the
resulting ``.aln`` file in at least two places (``protein_type`` /
interface-composition analysis, and ``swiss_model``) with subtly different,
copy-pasted parsing code. This module provides one alignment runner and one
parser used by both.
"""

from __future__ import annotations

import os
from dataclasses import dataclass
from pathlib import Path

from Bio.Align.Applications import ClustalwCommandline

from . import config


@dataclass
class PairwiseAlignment:
    """A gapped pairwise alignment between a monomer query and a mer candidate."""

    query_aligned: str
    mer_aligned: str

    def residue_index_map(self) -> tuple[dict[int, str], dict[int, str]]:
        """
        Build 1-based residue-number -> one-letter-code maps for both
        sequences, skipping alignment columns where the mer sequence has a
        gap (mirrors the original script's ``mono_dic`` / ``mer_dic``
        construction, which numbers residues by their position in the
        *mer* sequence).
        """
        query_map: dict[int, str] = {}
        mer_map: dict[int, str] = {}
        position = 1
        for query_char, mer_char in zip(self.query_aligned, self.mer_aligned):
            if mer_char == "-":
                continue
            mer_map[position] = mer_char
            query_map[position] = query_char
            position += 1
        return query_map, mer_map


def write_pair_fasta(query_fasta: Path, mer_fasta: Path, out_path: Path) -> Path:
    """Concatenate two 2-line FASTA files into a single multi-FASTA input for ClustalW."""
    query_lines = query_fasta.read_text().split("\n")
    mer_lines = mer_fasta.read_text().split("\n")
    out_path.parent.mkdir(parents=True, exist_ok=True)
    out_path.write_text(
        "\n".join([query_lines[0], query_lines[1], "", mer_lines[0], mer_lines[1], ""])
    )
    return out_path


def run_clustalw(fasta_path: Path) -> Path:
    """Run ClustalW on a multi-FASTA file and return the resulting ``.aln`` path."""
    command = ClustalwCommandline(config.CLUSTALW_EXE, infile=str(fasta_path))
    os.system(str(command))
    return fasta_path.with_suffix(".aln")


def parse_pairwise_aln(aln_path: Path) -> PairwiseAlignment:
    """
    Parse a 2-sequence ClustalW ``.aln`` file into a :class:`PairwiseAlignment`.

    ClustalW's block format interleaves the two sequences with a consensus
    line and a blank line; the original script picked out lines by
    ``count % 4`` -- kept here, but named.
    """
    lines = aln_path.read_text().split("\n")[3:]
    query_seq, mer_seq = "", ""
    for i, line in enumerate(lines):
        parts = line.split()
        if not parts or len(parts) < 2:
            continue
        if i % 4 == 0:
            query_seq += parts[1]
        elif i % 4 == 1:
            mer_seq += parts[1]
    return PairwiseAlignment(query_aligned=query_seq, mer_aligned=mer_seq)


def align_pair(query_fasta: Path, mer_fasta: Path, workdir: Path, label: str) -> PairwiseAlignment:
    """High-level helper: write a pair FASTA, run ClustalW, and parse the result."""
    pair_fasta = write_pair_fasta(query_fasta, mer_fasta, workdir / f"{label}.fasta")
    aln_path = run_clustalw(pair_fasta)
    return parse_pairwise_aln(aln_path)
