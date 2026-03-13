"""Pytest helpers for keeping generated test artifacts out of the repo tree."""

from pathlib import Path

import pytest


REPO = Path(__file__).resolve().parents[1]
ARTIFACT_DIRS = (REPO, REPO / "tests")
ARTIFACT_SUFFIXES = {
    ".dot",
    ".gml",
    ".jpg",
    ".log",
    ".pdf",
    ".png",
    ".py",
    ".sh",
    ".svg",
    ".tex",
    ".tsv",
    ".txt",
    ".xlsx",
}


def _iter_artifacts():
    """Yield top-level generated files commonly produced by the test suite."""
    for directory in ARTIFACT_DIRS:
        if not directory.exists():
            continue
        for path in directory.iterdir():
            if path.is_file() and path.suffix.lower() in ARTIFACT_SUFFIXES:
                yield path


@pytest.fixture(autouse=True, scope="session")
def cleanup_generated_artifacts():
    """Remove new top-level artifacts created by the test session."""
    before = set(_iter_artifacts())
    yield
    for path in set(_iter_artifacts()).difference(before):
        path.unlink(missing_ok=True)
