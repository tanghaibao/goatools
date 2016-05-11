#!/usr/bin/env python
"""This test addresses "Error printing GOTerm #61"."""

import sys
import goatools

def test_go_print(prt=sys.stdout):
    """Test that all GO Terms can be printed, even if level/depth are not assigned."""
    reader = goatools.obo_parser.OBOReader('go-basic.obo')
    prt.write("\n{OBJ}\n\n".format(OBJ=reader))
    go_terms = list(reader)
    prt.write("First GO Record: {REC}\n".format(REC=go_terms[0]))
    for idx, go_rec in enumerate([rec for rec in reader]):
        prt.write("{I:>7,} {RECORD}\n".format(I=idx, RECORD=go_rec))

if __name__ == '__main__':
  test_go_print()
