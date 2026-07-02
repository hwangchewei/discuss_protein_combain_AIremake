"""
Selenium WebDriver helpers.

Several sites used by this pipeline (NCBI VAST+, EBI PISA, the DALI server,
SWISS-MODEL) either require JavaScript interaction or don't expose a clean
API, so the original script drove a real Chrome browser via Selenium. This
module centralises driver creation (the original repeated
``webdriver.Chrome(executable_path='./chromedriver', chrome_options=...)``
-- a pattern removed in modern Selenium -- in every function that needed a
browser) and a couple of small polling helpers.
"""

from __future__ import annotations

import contextlib
import logging
from typing import Iterator

from selenium import webdriver
from selenium.webdriver.chrome.options import Options
from selenium.webdriver.chrome.service import Service

from . import config

logger = logging.getLogger(__name__)


def build_chrome_options(*, headless: bool = True, download_dir: str | None = None) -> Options:
    """Build a Chrome ``Options`` object with the pipeline's usual settings."""
    options = Options()
    if headless:
        options.add_argument("--headless=new")
        options.add_argument("--disable-gpu")
    options.add_argument("log-level=3")
    options.add_argument(f"user-agent={config.USER_AGENT_STRING}")
    if download_dir:
        options.add_experimental_option("prefs", {"download.default_directory": download_dir})
    return options


def new_driver(*, headless: bool = True, download_dir: str | None = None) -> webdriver.Chrome:
    """
    Create a new Chrome WebDriver instance.

    Uses ``selenium.webdriver.chrome.service.Service`` instead of the
    deprecated ``executable_path=`` keyword argument used throughout the
    original script.
    """
    service = Service(executable_path=config.CHROMEDRIVER_PATH)
    options = build_chrome_options(headless=headless, download_dir=download_dir)
    return webdriver.Chrome(service=service, options=options)


@contextlib.contextmanager
def chrome_driver(*, headless: bool = True, download_dir: str | None = None) -> Iterator[webdriver.Chrome]:
    """Context manager that guarantees ``driver.quit()`` is called on exit."""
    driver = new_driver(headless=headless, download_dir=download_dir)
    try:
        yield driver
    finally:
        driver.quit()
