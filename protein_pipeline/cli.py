"""
Command-line entry point.

The original ``main()`` was a single function that interleaved three
different, non-reusable workflows and prompted interactively for paths and
thresholds with bare ``input()`` calls -- undocumented, unscriptable, and
impossible to re-run partially. This module exposes each stage as its own
subcommand so the pipeline can be run, resumed, or scripted stage by stage.

Example
-------
    python -m protein_pipeline.cli download-vast monomer.txt
    python -m protein_pipeline.cli filter-neighbours 1ABC monomer.txt homomer.txt
    python -m protein_pipeline.cli interface-report
"""

from __future__ import annotations

import argparse
import logging
import sys
from multiprocessing import Pool
from pathlib import Path

from . import (
    config,
    dali,
    fasta_utils,
    interface_composition,
    pisa,
    reporting,
    sorting,
    swissmodel,
    vast,
)

logger = logging.getLogger(__name__)


def _read_pdb_id_column(path: Path) -> list[str]:
    """Read PDB ids from column 2 of a whitespace-delimited list file."""
    ids = []
    for line in path.read_text().splitlines():
        parts = line.split()
        if len(parts) > 1:
            ids.append(parts[1].upper())
    return ids


def cmd_download_vast(args: argparse.Namespace) -> None:
    """Stage 1: fetch VAST+ JSON for every monomer id, in parallel."""
    pdb_ids = _read_pdb_id_column(Path(args.monomer_list))
    with Pool(args.workers) as pool:
        pool.map(vast.download_vastplus_json_via_browser_or_http, pdb_ids)


def cmd_filter_neighbours(args: argparse.Namespace) -> None:
    """Stage 2: filter VAST+ neighbours into candidate homomer FASTAs."""
    vast.collect_candidate_fastas(
        args.pdb_id, Path(args.monomer_list), Path(args.homomer_list)
    )
    vast.sort_by_oligomeric_state(args.pdb_id)
    fasta_utils.fetch_missing_fastas(args.pdb_id)
    fasta_utils.filter_aligned_neighbours(args.pdb_id)


def cmd_pisa_stage(args: argparse.Namespace) -> None:
    """Stage 3: download & filter PISA interface data, strip ligands."""
    pisa.download_interfaces_xml(args.pdb_id)
    pisa.filter_true_protein_interfaces(args.pdb_id)
    pisa.check_already_known_assembly(args.pdb_id)
    pisa.move_promoted_files(args.pdb_id)
    pisa.strip_ligands(args.pdb_id)
    pisa.check_already_known_assembly(args.pdb_id)


def cmd_dali_stage(args: argparse.Namespace) -> None:
    """Stage 4: reconstruct candidate assemblies via the DALI server."""
    dali.process_all_candidates(args.pdb_id)


def cmd_swissmodel_stage(args: argparse.Namespace) -> None:
    """Stage 5: build homology models for remaining candidates."""
    swissmodel.build_models_for_all_candidates(args.pdb_id)


def cmd_interface_report(args: argparse.Namespace) -> None:
    """Stage 6: compute per-interface amino-acid composition."""
    interface_composition.whole_sequence_composition(args.pdb_id)
    interface_composition.interface_residue_composition(args.pdb_id)


def cmd_aggregate_reports(args: argparse.Namespace) -> None:
    """Stage 7: build the aggregate change/composition CSV reports."""
    pdb_ids = reporting.interface_change_report()
    reporting.interface_amino_change_report()
    print(f"{len(pdb_ids)} PDB entries with recorded interfaces")


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="protein_pipeline",
        description="Homomer-interface discovery and composition analysis pipeline.",
    )
    parser.add_argument(
        "--root", type=Path, default=None,
        help="Override the pipeline data root (default: %s or PROTEIN_PIPELINE_ROOT env var)" % config.PROJECT_ROOT,
    )
    subparsers = parser.add_subparsers(dest="command", required=True)

    p = subparsers.add_parser("download-vast", help="Fetch VAST+ JSON for a list of monomer PDB ids")
    p.add_argument("monomer_list", help="Path to a monomer id list file (PDBePISA format)")
    p.add_argument("--workers", type=int, default=6)
    p.set_defaults(func=cmd_download_vast)

    p = subparsers.add_parser("filter-neighbours", help="Filter VAST+ neighbours into candidate FASTAs")
    p.add_argument("pdb_id")
    p.add_argument("monomer_list")
    p.add_argument("homomer_list")
    p.set_defaults(func=cmd_filter_neighbours)

    p = subparsers.add_parser("pisa-stage", help="Download & filter PISA interface data")
    p.add_argument("pdb_id")
    p.set_defaults(func=cmd_pisa_stage)

    p = subparsers.add_parser("dali-stage", help="Reconstruct candidate assemblies via DALI")
    p.add_argument("pdb_id")
    p.set_defaults(func=cmd_dali_stage)

    p = subparsers.add_parser("swissmodel-stage", help="Build homology models for candidates")
    p.add_argument("pdb_id")
    p.set_defaults(func=cmd_swissmodel_stage)

    p = subparsers.add_parser("interface-report", help="Compute per-interface amino-acid composition")
    p.add_argument("pdb_id")
    p.set_defaults(func=cmd_interface_report)

    p = subparsers.add_parser("aggregate-reports", help="Build aggregate change/composition CSV reports")
    p.set_defaults(func=cmd_aggregate_reports)

    return parser


def main(argv: list[str] | None = None) -> int:
    logging.basicConfig(level=logging.INFO, format="%(levelname)s %(name)s: %(message)s")
    parser = build_parser()
    args = parser.parse_args(argv)

    if args.root:
        config.DATA_DIR = args.root  # type: ignore[misc]

    args.func(args)
    return 0


if __name__ == "__main__":
    sys.exit(main())
