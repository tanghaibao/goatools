#!/usr/bin/env python
"""Test CompareGOsCli"""

import os
from goatools.cli.compare_gos import CompareGOsCli

__copyright__ = "Copyright (C) 2016-2019, DV Klopfenstein, H Tang, All rights reserved."

REPO = os.path.join(os.path.dirname(os.path.abspath(__file__)), "../")


def test_example():
    """Test CompareGOsCli"""
    files = ['tat_gos_simple1.tsv', 'tat_gos_simple2.tsv']
    kws = {
        'GO_FILE': [os.path.join(REPO, 'data/compare_gos', f) for f in files],
        'obo': os.path.join(REPO, 'go-basic.obo'),
        'slims': os.path.join(REPO, 'goslim_generic.obo'),
        'sections': os.path.join(REPO, 'data/compare_gos/sections.txt'),
    }
    obj = CompareGOsCli(**kws)
    obj.write(fout_xlsx='compare_gos_v1.xlsx', fout_txt='compare_gos_v1.txt', verbose=True)
    obj.write(fout_xlsx='compare_gos_v0.xlsx', fout_txt='compare_gos_v0.txt', verbose=False)


if __name__ == '__main__':
    test_example()

# Copyright (C) 2016-2019, DV Klopfenstein, H Tang, All rights reserved.
