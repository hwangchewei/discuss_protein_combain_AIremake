# protein_pipeline

A modular, re-engineered rebuild of the original `discuss_protein_combain.py`
/ `__init__.py` scripts: a research pipeline that starts from a list of
monomeric PDB entries, finds homomeric structural neighbours (NCBI VAST+),
downloads and aligns their sequences/structures (RCSB, EBI PISA, the DALI
server, SWISS-MODEL), and characterises the resulting protein-protein
interfaces by amino-acid composition (hydrophobic / hydrophilic / charged).

## What changed from the original scripts

The original code was a single ~1,700-line file plus a ~500-line helper
module, both hard-coded to Windows drive paths (`D:\...`, `E:\...`), full of
duplicated logic, 20-branch `if/elif` chains for amino-acid tallying, dead
commented-out code, bare `except:` clauses, and `input()`-driven interactive
prompts baked into `main()`. This rebuild keeps the same overall research
workflow but:

- splits it into 14 focused modules (see below) instead of 2 monolithic files;
- replaces hard-coded paths with a single configurable `config.py`;
- replaces `executable_path=` (deprecated in modern Selenium) with
  `selenium.webdriver.chrome.service.Service`;
- consolidates 4+ separately-defined copies of the amino-acid
  classification tables into one `amino_acids.py`;
- replaces 20-branch `if/elif` letter-counting with `collections.Counter`;
- removes dead/commented-out code blocks (most notably in the surface
  analysis class, which had several superseded scoring formulas commented
  out but left in place);
- exposes each pipeline stage as its own CLI subcommand instead of one
  monolithic, interactively-prompted `main()`.

## Module layout

| Module | Responsibility | Replaces |
|---|---|---|
| `config.py` | paths, thresholds, network/tool settings | scattered hard-coded paths & constants |
| `amino_acids.py` | amino-acid codes & physicochemical classification | duplicated tables in both original files |
| `http_client.py` | plain HTTP downloads (FASTA, VAST+ JSON, PISA XML, PDB) | repeated `requests.get(...).text` blocks |
| `browser.py` | Selenium WebDriver setup | repeated `webdriver.Chrome(...)` blocks |
| `structure_surface.py` | Bio.PDB surface/interface geometry (`SurfaceAnalyzer`) | `__init__.py`'s `biopython` class |
| `vast.py` | NCBI VAST+ neighbour discovery & filtering | `json_download`, `no_selection`, `monomer_fasta`, `oligState` |
| `fasta_utils.py` | FASTA retrieval & alignment-ratio filtering | `get_fasta`, `alignedResidues`, `not_monomer` |
| `alignment.py` | ClustalW pairwise alignment + residue-index mapping | inline ClustalW calls in `protein_type`/`swiss_model` |
| `pisa.py` | EBI PISA interface downloads & interactive checks | `pisa_download`, `check_mer`, `pisa_mer_check(_)`, `pdb_del_ligand`, `move_file`, `pisa` |
| `dali.py` | DALI server structural superposition | `dali`, `dali_start` |
| `swissmodel.py` | SWISS-MODEL homology modelling | `swiss_model` |
| `interface_composition.py` | per-interface / per-protein amino-acid composition | `protein_type`, `interface_check` |
| `reporting.py` | aggregate CSV/text reports across many interfaces | `interface_change`, `interface_amino_change`, `protein_face` |
| `sorting.py` | misc. bucket-sorting housekeeping | `is_dimer`, `check_fasta_pdb`, `pdb_not_in_file` |
| `cli.py` | command-line entry point | `main()` |

## Requirements

See `requirements.txt`. You will also need, on `PATH` (or pointed to via the
environment variables below):

- Google Chrome + a matching `chromedriver` (for the VAST+/PISA/DALI/
  SWISS-MODEL stages, which drive a real browser)
- `clustalw2` (for pairwise sequence alignment)

## Configuration

All paths derive from a single root, configurable via an environment
variable instead of hard-coded drive letters:

```bash
export PROTEIN_PIPELINE_ROOT=/path/to/monomer_json
export CHROMEDRIVER_PATH=/path/to/chromedriver   # optional, default: "chromedriver" on PATH
export CLUSTALW_EXE=/path/to/clustalw2           # optional, default: "clustalw2" on PATH
```

## Usage

```bash
# Stage 1: fetch VAST+ neighbour data for a list of monomer PDB ids
python -m protein_pipeline.cli download-vast monomer.txt

# Stage 2: filter neighbours into candidate homomer FASTAs for one entry
python -m protein_pipeline.cli filter-neighbours 1ABC monomer.txt homomer.txt

# Stage 3: PISA interface download & filtering
python -m protein_pipeline.cli pisa-stage 1ABC

# Stage 4/5: structural reconstruction / homology modelling
python -m protein_pipeline.cli dali-stage 1ABC
python -m protein_pipeline.cli swissmodel-stage 1ABC

# Stage 6: per-interface amino-acid composition
python -m protein_pipeline.cli interface-report 1ABC

# Stage 7: aggregate reports across everything analysed so far
python -m protein_pipeline.cli aggregate-reports
```

Each stage can also be called directly as a Python function -- see the
module docstrings for the mapping back to the original script's functions.

## Notes / known limitations carried over from the original design

- The PISA/DALI/SWISS-MODEL automation drives real websites via Selenium
  and is inherently coupled to their current page structure (long XPath
  expressions). If those sites change their HTML, the corresponding
  functions in `pisa.py`, `dali.py`, and `swissmodel.py` will need updating.
- This pipeline was not runnable in the sandbox used to produce this
  rebuild (no network access, no Selenium/Biopython installed), so the
  network- and browser-driving code paths are refactored for clarity and
  correctness of *intent* but have not been exercised end-to-end. Please
  test each stage against your own environment before relying on it.
