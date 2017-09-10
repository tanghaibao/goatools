#!/usr/bin/env python
"""This test addresses "Error printing GOTerm #61"."""

import os
import sys
import goatools
from goatools.base import download_go_basic_obo

def test_go_print(prt=sys.stdout):
    """Test that all GO Terms can be printed, even if level/depth are not assigned."""
    prt_pypath(prt)
    file_obo = os.path.join(os.getcwd(), "go-basic.obo")
    obo_file = download_go_basic_obo(file_obo, prt=prt, loading_bar=None)
    reader = goatools.obo_parser.OBOReader(obo_file)
    go_terms = list(reader)
    prt.write("Python Version: {VER}\n\n".format(VER=sys.version))
    prt.write("\nOBOReader: {OBJ}\n\n".format(OBJ=reader))
    prt.write("format-version: {VER}\n".format(VER=reader.format_version))
    prt.write("data-version: {VER}\n\n".format(VER=reader.data_version))
    prt.write("Found {N} GO Records:\n".format(N=len(go_terms)))
    for idx, go_rec in enumerate(go_terms):
        prt.write("{I:>7,} {RECORD}\n".format(I=idx, RECORD=go_rec))

def prt_pypath(prt):
    """Print PYTHONPATH contents."""
    pypathes = os.environ.get('PYTHONPATH', None)
    if pypathes:
        prt.write("\nPYTHONPATH:\n")
        for idx, pypath in enumerate(pypathes.split(os.pathsep)):
            prt.write("    {IDX} {PATH}\n".format(IDX=idx, PATH=pypath))
        prt.write("\n")


if __name__ == '__main__':
    test_go_print()
