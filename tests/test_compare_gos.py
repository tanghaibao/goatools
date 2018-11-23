#!/usr/bin/env python
"""Test gracefully exiting if no study genes are in assc or population."""

import os
from goatools.cli.compare_gos import CompareGOs

__copyright__ = "Copyright (C) 2016-2019, DV Klopfenstein, H Tang, All rights reserved."

REPO = os.path.join(os.path.dirname(os.path.abspath(__file__)), "../")


def test_example():
    """Test GOEnrichmentStudy::print_results."""
    files = ['tat_gos_simple1.tsv', 'tat_gos_simple2.tsv']
    kws = {
        'GO_FILE': [os.path.join(REPO, 'data/compare_gos', f) for f in files],
        'obo': os.path.join(REPO, 'go-basic.obo'),
        'slims': os.path.join('goslim_generic.obo'),
        'sections': os.path.join('data/compare_gos/sections.txt'),
        'ofile': 'compare_gos.txt',
        'xlsx': 'compare_gos.xlsx',
    }
    obj = CompareGOs(**kws)
    obj.cli()


if __name__ == '__main__':
    test_example()

# Copyright (C) 2016-2019, DV Klopfenstein, H Tang, All rights reserved.
