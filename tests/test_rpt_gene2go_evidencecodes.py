#!/usr/bin/env python3
"""Tests that all evidence codes seen in NCBI's gene2go have description."""

from __future__ import print_function

__copyright__ = "Copyright (C) 2016-2019, DV Klopfenstein, H Tang. All rights reserved."
__author__ = "DV Klopfenstein"

import os

from goatools.associations import dnld_ncbi_gene_file
from goatools.evidence_codes import EvidenceCodes
REPO = os.path.join(os.path.dirname(os.path.abspath(__file__)), "../")

def test_ev():
    """Return GO associations from a GAF file. Download if necessary."""
    evs = _get_evidencecodes('gene2go')
    obj = EvidenceCodes()
    missing = evs.difference(obj.code2nt)
    assert not missing, 'MISSING({EV})'.format(EV=missing)

def _get_evidencecodes(fin_gene2go):
    """Get all evidence codes and qualifiers."""
    evs = set()
    fin_gene2go = os.path.join(REPO, 'gene2go')
    dnld_ncbi_gene_file(fin_gene2go, force_dnld=False, loading_bar=False)
    with open(fin_gene2go) as ifstrm:
        for line in ifstrm:
            if line[0] != '#': # Line contains data. Not a comment
                line = line.rstrip() # chomp
                flds = line.split('\t')
                if len(flds) >= 5:
                    # taxid_curr, geneid, go_id, evidence, qualifier = flds[:5]
                    evidence = flds[3]
                    assert len(evidence) >= 2, flds
                    evs.add(evidence)
    print('{N} evidence codes in {FIN}'.format(N=len(evs), FIN=fin_gene2go))
    return evs


if __name__ == '__main__':
    test_ev()

# Copyright (C) 2016-2019, DV Klopfenstein, H Tang. All rights reserved.
