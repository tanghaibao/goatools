#!/usr/bin/env python3
"""Test commands run in notebooks/get_isa_and_partof.ipynb"""

from __future__ import print_function

__copyright__ = "Copyright (C) 2010-2019, DV Klopfenstein, H Tang. All rights reserved."

import os
import sys

# REPO = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..")


def test_plotgos(run_all=False):
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
        'python3 scripts/go_plot.py -o viral_r0.png                                                GO:0019222#d8dcd6 GO:0060150 --obo=tests/data/i126/viral_gene_silence.obo --go_color_file=tests/data/i126/viral_gene_silence.txt',
        'python3 scripts/go_plot.py -o viral_r1.png -r                                             GO:0010468#d8dcd6 GO:0060150 --obo=tests/data/i126/viral_gene_silence.obo --go_color_file=tests/data/i126/viral_gene_silence.txt',
        'python3 scripts/go_plot.py -o viral_r_partof.png --relationships=part_of                  GO:0010468#d8dcd6 GO:0060150 --obo=tests/data/i126/viral_gene_silence.obo --go_color_file=tests/data/i126/viral_gene_silence.txt',
        'python3 scripts/go_plot.py -o viral_reg.png --relationships=regulates                     GO:0050794#d8dcd6 GO:0060150 --obo=tests/data/i126/viral_gene_silence.obo --go_color_file=tests/data/i126/viral_gene_silence.txt',
        'python3 scripts/go_plot.py -o viral_rp.png --relationships=positively_regulates           GO:0048522#d8dcd6 GO:0060150 --obo=tests/data/i126/viral_gene_silence.obo --go_color_file=tests/data/i126/viral_gene_silence.txt',
        'python3 scripts/go_plot.py -o viral_rn.png --relationships=regulates,negatively_regulates GO:0050794#d8dcd6 GO:0060150 --obo=tests/data/i126/viral_gene_silence.obo --go_color_file=tests/data/i126/viral_gene_silence.txt',
    ]


if __name__ == '__main__':
    test_plotgos(len(sys.argv) != 1)

# Copyright (C) 2010-2019, DV Klopfenstein, H Tang. All rights reserved.
