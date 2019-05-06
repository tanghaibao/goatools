#!/usr/bin/env python
"""Tests read a GAF with missing (required) DB_Symbol text."""

import os
from goatools.anno.gaf_reader import GafReader

REPO = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..")


def test_missingsym():
    """Tests read a GAF with missing (required) DB_Symbol text."""
    # Original gaf file (mgi.gaf) was reduced
    fin_gaf = "tests/data/gaf_missingsym.mgi"
    # Test that gene products that are missing the required DB_Symbol are ignored
    gafobj = GafReader(fin_gaf, hdr_only=False)
    assert not gafobj.chk_associations('gaf_missingsym.err')

    # assert len(gene2gos) == 16, len(gene2gos)
    # assert 'MGI:3643263' not in gene2gos
    # assert 'P84751' not in gene2gos
    # # Test that gene products that are missing the required DB_Symbol are ignored
    # gene2gos = read_gaf(os.path.join(REPO, fin_gaf), allow_missing_symbol=False)
    # assert len(gene2gos) == 16, len(gene2gos)
    # assert 'MGI:3643263' not in gene2gos
    # assert 'P84751' not in gene2gos
    # # Tests saving annotation, even if missing required DB_Symbol
    # gene2gos = read_gaf(os.path.join(REPO, fin_gaf), allow_missing_symbol=True)
    # assert len(gene2gos) == 18
    # assert 'MGI:3643263' in gene2gos
    # assert 'P84751' in gene2gos


if __name__ == '__main__':
    test_missingsym()
