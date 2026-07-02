"""
Plain (non-Selenium) HTTP downloads.

The original script re-implemented ``requests.get(url=..., headers=...).text``
followed by writing to disk in half a dozen places (``no_selection``,
``monomer_fasta``, ``get_fasta``, ``pisa_download``, ``dali_start``...) with
no error handling. This module centralises those calls with retries and
clear function names.
"""

from __future__ import annotations

import logging
from pathlib import Path

import requests
from requests.adapters import HTTPAdapter
from urllib3.util.retry import Retry

from . import config

logger = logging.getLogger(__name__)

_session: requests.Session | None = None


def get_session() -> requests.Session:
    """Return a shared :class:`requests.Session` with sane retry behaviour."""
    global _session
    if _session is None:
        _session = requests.Session()
        _session.headers.update(config.REQUEST_HEADERS)
        retries = Retry(total=3, backoff_factor=1.0, status_forcelist=(429, 500, 502, 503, 504))
        _session.mount("https://", HTTPAdapter(max_retries=retries))
        _session.mount("http://", HTTPAdapter(max_retries=retries))
    return _session


def fetch_text(url: str) -> str:
    """GET a URL and return the response body as text, raising on HTTP errors."""
    response = get_session().get(url, timeout=config.REQUEST_TIMEOUT)
    response.raise_for_status()
    return response.text


def fetch_vastplus_json(pdb_id: str) -> str:
    """Fetch the raw NCBI VAST+ JSON payload for a PDB id."""
    url = f"https://www.ncbi.nlm.nih.gov/Structure/vastplus/vastplus.cgi?uid={pdb_id}&getdata=json"
    return fetch_text(url)


def fetch_rcsb_fasta(pdb_id: str) -> str:
    """Fetch the FASTA sequence(s) for a PDB entry from RCSB."""
    url = f"https://www.rcsb.org/fasta/entry/{pdb_id}/display"
    return fetch_text(url)


def fetch_pisa_interfaces_xml(pdb_id: str) -> str:
    """Fetch the PISA interface-analysis XML for a PDB entry."""
    url = f"https://www.ebi.ac.uk/pdbe/pisa/cgi-bin/interfaces.pisa?{pdb_id.lower()}"
    return fetch_text(url)


def fetch_pdb_file(pdb_id: str) -> str:
    """Fetch the raw PDB-format coordinate file for a PDB entry."""
    url = f"https://files.rcsb.org/view/{pdb_id}.pdb"
    return fetch_text(url)


def save_text(path: Path, content: str) -> None:
    """Write text content to ``path``, creating parent directories as needed."""
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(content, encoding="utf-8")


def download_to_file(fetch_fn, identifier: str, destination: Path) -> Path:
    """
    Fetch content via ``fetch_fn(identifier)`` and write it to ``destination``.

    Small helper that replaces the copy-pasted
    ``fasta = requests.get(...).text; open(path, 'w').write(fasta)`` pattern.
    """
    content = fetch_fn(identifier)
    save_text(destination, content)
    return destination
