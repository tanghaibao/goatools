#!/usr/bin/env python
"""Test validity checks for the id2gos annotation format (reader_idtogos.py)."""
# pylint: disable=line-too-long

import os
import tempfile

import pytest

from goatools.anno.idtogos_reader import IdToGosReader

__copyright__ = (
    "Copyright (C) 2016-present, DV Klopfenstein, H Tang. All rights reserved."
)

REPO = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..")

# ---------------------------------------------------------------------------
# Helper: write a temporary id2gos file and return its path
# ---------------------------------------------------------------------------

def _write_tmp(lines):
    """Write *lines* to a temp file and return the file path."""
    fout = tempfile.NamedTemporaryFile(
        mode="w", suffix=".txt", delete=False
    )
    fout.write("\n".join(lines) + "\n")
    fout.close()
    return fout.name


# ---------------------------------------------------------------------------
# Fixtures / shared data
# ---------------------------------------------------------------------------

VALID_LINES = [
    "GENE1\tGO:0005575;GO:0003674;GO:0006970",
    "GENE2\tGO:0005575;GO:0003674;GO:0040029",
]

# Lines that use "; " (semicolon+space) instead of ";" – the bug in the issue
SPACE_SEP_LINES = [
    "GENE1\tGO:0005575; GO:0003674; GO:0006970",
    "GENE2\tGO:0005575; GO:0003674; GO:0040029",
]

# Lines that contain a deliberately malformed GO ID
MALFORMED_LINES = [
    "GENE1\tGO:0005575;NOTAGO:123;GO:0006970",
    "GENE2\tGO:0005575;GO:0003674",
]


# ---------------------------------------------------------------------------
# Tests
# ---------------------------------------------------------------------------

def test_valid_file_no_warnings():
    """A valid id2gos file should be read with all GO IDs loaded."""
    fin = _write_tmp(VALID_LINES)
    try:
        obj = IdToGosReader(fin)
        id2gos = obj.get_id2gos_nss()
        assert len(id2gos) == 2, id2gos
        # chk_associations should find no errors
        assert obj.chk_associations() == 0
    finally:
        os.unlink(fin)


def test_space_separated_goids_are_stripped():
    """GO IDs with surrounding spaces (e.g. from '; ' separator) are stripped
    and accepted, not silently discarded."""
    fin = _write_tmp(SPACE_SEP_LINES)
    try:
        obj = IdToGosReader(fin)
        id2gos = obj.get_id2gos_nss()
        # Both genes must be present
        assert len(id2gos) == 2, id2gos
        # Each gene must have all 3 GO IDs (spaces stripped before validation)
        for gene_id, gos in id2gos.items():
            assert len(gos) == 3, (
                f"Expected 3 GO IDs for {gene_id}, got {len(gos)}: {gos}"
            )
        # chk_associations should find no errors (file has spaces stripped on read)
        assert obj.chk_associations() == 0
    finally:
        os.unlink(fin)


def test_malformed_goid_is_excluded():
    """A malformed GO ID should be excluded from the loaded associations."""
    fin = _write_tmp(MALFORMED_LINES)
    try:
        obj = IdToGosReader(fin)
        id2gos = obj.get_id2gos_nss()
        # GENE1 has one bad GO ID → only 2 valid IDs should survive
        assert "GO:0005575" in id2gos.get("GENE1", set()), id2gos
        assert "NOTAGO:123" not in id2gos.get("GENE1", set()), id2gos
        assert len(id2gos["GENE1"]) == 2, id2gos["GENE1"]
        # chk_associations should report exactly 1 invalid ID
        assert obj.chk_associations() == 1
    finally:
        os.unlink(fin)


def test_chk_associations_writes_error_file():
    """chk_associations should write invalid IDs to the error file when requested."""
    fin = _write_tmp(MALFORMED_LINES)
    ferr = tempfile.NamedTemporaryFile(
        mode="w", suffix=".err", delete=False
    ).name
    try:
        obj = IdToGosReader(fin)
        num_invalid = obj.chk_associations(fout_err=ferr)
        assert num_invalid == 1
        with open(ferr) as fh:
            content = fh.read()
        assert "NOTAGO:123" in content
        assert "LINE" in content
    finally:
        os.unlink(fin)
        if os.path.exists(ferr):
            os.unlink(ferr)


def test_small_association_no_errors():
    """The bundled small_association test file should contain no invalid GO IDs."""
    fin = os.path.join(REPO, "tests", "data", "small_association")
    assert os.path.exists(fin), f"Test data not found: {fin}"
    obj = IdToGosReader(fin)
    assert obj.chk_associations() == 0


if __name__ == "__main__":
    pytest.main([__file__, "-v"])
