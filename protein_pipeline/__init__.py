"""
protein_pipeline
=================

A modular, re-engineered version of the original monolithic
``discuss_protein_combain.py`` script.

The pipeline automates a research workflow that, starting from a list of
monomeric PDB entries, looks for homomeric neighbours (via NCBI VAST+),
downloads sequences/structures (RCSB, PISA, DALI, SWISS-MODEL), aligns
sequences (ClustalW), and characterises the resulting protein-protein
interfaces in terms of amino-acid composition (hydrophobic / hydrophilic /
positively or negatively charged).

The original script mixed all of this into a single 1,700-line file with
hard-coded Windows paths, global mutable state, duplicated logic, and no
error handling. This package keeps the same overall workflow but splits it
into focused modules:

- ``config``               : central configuration (paths, constants, user agent)
- ``amino_acids``          : amino-acid codes / classification helpers
- ``http_client``          : plain HTTP downloads (FASTA, VAST+ JSON, PISA XML, PDB)
- ``browser``              : Selenium WebDriver helpers
- ``structure_surface``    : Bio.PDB based surface/interface geometry analysis
- ``vast``                 : NCBI VAST+ neighbour discovery & filtering
- ``fasta_utils``          : FASTA retrieval & alignment-based sequence filtering
- ``alignment``            : ClustalW alignment + monomer<->mer residue mapping
- ``pisa``                 : EBI PISA interface downloads & interactive checks
- ``dali``                 : DALI server structural superposition automation
- ``swissmodel``           : SWISS-MODEL homology modelling automation
- ``interface_composition``: per-interface amino-acid composition analysis
- ``reporting``            : aggregate CSV/text report generation
- ``cli``                  : command line entry point tying the stages together

Only orchestration (``cli.py``) and truly generic helpers are imported here
to keep import time low and avoid requiring optional dependencies (Selenium,
Biopython) unless the relevant stage is actually used.
"""

__version__ = "1.0.0"

__all__ = ["__version__"]
