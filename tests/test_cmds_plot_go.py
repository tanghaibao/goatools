#!/usr/bin/env python3
"""Test running an enrichment using any annotation file format."""

from __future__ import print_function

__copyright__ = "Copyright (C) 2010-2019, DV Klopfenstein, H Tang. All rights reserved."

import os
import sys

# REPO = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..")


def test_cmds_plot_annos(run_all=False):
    """RUn an enrichments using all annotation file formats"""

    if run_all:
        for idx, cmd in enumerate(_get_cmds()):
            print('------------------- TEST {I} ------------------------------------'.format(I=idx))
            print('CMD: {CMD}'.format(CMD=cmd))
            assert os.system(cmd) == 0
        print("TEST PASSED")
    else:
        print('RUN THIS TEST WITH AN ARGUMENT')


def _get_cmds():
    """Get commands used in ./doc/md/README_find_enrichment.md"""
    # pylint: disable=line-too-long
    return [
        'scripts/go_plot.py GO:0016150 -o aaaa_none.png',
        'scripts/go_plot.py GO:0016150 -o aaaa_gpad.png --gpad goa_human.gpa',
        'scripts/go_plot.py GO:0016150 -o aaaa_gaf.png --gaf goa_human.gaf',
        'scripts/go_plot.py GO:0016150 -o aaaa_orig.png --id2gos data/association',
    ]


if __name__ == '__main__':
    test_cmds_plot_annos(len(sys.argv) != 1)

# Copyright (C) 2010-2019, DV Klopfenstein, H Tang. All rights reserved.
