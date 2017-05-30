#!/usr/bin/env python
"""This test addresses "Error printing GOTerm #61"."""

import sys
import goatools
from goatools.base import download_go_basic_obo

def test_go_print(prt=sys.stdout):
    """Test that all GO Terms can be printed, even if level/depth are not assigned."""
    obo_file = download_go_basic_obo(prt=prt)
    reader = goatools.obo_parser.OBOReader(obo_file)
    go_terms = list(reader)
    prt.write("\nOBOReader: {OBJ}\n\n".format(OBJ=reader))
    prt.write("format-version: {VER}\n".format(VER=reader.format_version))
    prt.write("data-version: {VER}\n\n".format(VER=reader.data_version))
    prt.write("Found {N} GO Records:\n".format(N=len(go_terms)))
    for idx, go_rec in enumerate(go_terms):
        prt.write("{I:>7,} {RECORD}\n".format(I=idx, RECORD=go_rec))

if __name__ == '__main__':
    test_go_print()
