#!/usr/bin/env python
"""Test finding and reporting illegal GAF lines seen in the field."""

from __future__ import print_function

import os
import sys
from goatools.anno.gaf_reader import GafReader

REPO = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..")

def test_gaf_illegal(prt=sys.stdout):
    """Test finding and reporting illegal GAF lines seen in the field."""
    fin_gaf = os.path.join(REPO, 'data/gaf/goa_human_illegal.gaf')
    gafobj = GafReader(fin_gaf, hdr_only=False, prt=prt)
    # id2gos = gafobj.read_gaf()  # Read associations
    # for ntd in gafobj.associations:
    #     print(ntd)
    assert not gafobj.chk_associations('goa_human_illegal.err')
    #assert len(id2gos) == 15, "IDS FOUND: EXP({E}) ACT({A})".format(E=15, A=len(id2gos))
    #assert len(gafobj.datobj.ignored) == 1
    #assert len(gafobj.datobj.illegal_lines['ILLEGAL TAXON']) == 1


if __name__ == '__main__':
    test_gaf_illegal()
