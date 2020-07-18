#!/usr/bin/env python3
"""Create Table 2 (depth-01 child count) in GOATOOLS 2018 paper by DV Klopfenstein et al."""

__copyright__ = "Copyright (C) 2016-present, DV Klopfenstein, H Tang, All rights reserved."
__author__ = "DV Klopfenstein"

from goatools.base import get_godag
from goatools.gosubdag.godag_rcnt import CountRelatives
from goatools.gosubdag.rpt.wr_xlsx import GoDepth1LettersWr


def main():
    """Write Table of depth-01 GO terms w/child count"""

    fout_tex = "gos_depth01.tex"

    fin_dag = 'go-basic.obo'
    godag = get_godag(fin_dag, optional_attrs='relationship')
    rcntobj = CountRelatives(godag, relationships=True, dcnt=True)
    wrobj = GoDepth1LettersWr(rcntobj)
    wrobj.wr_tex(fout_tex)


if __name__ == '__main__':
    main()

# Copyright (C) 2016-present, DV Klopfenstein, H Tang, All rights reserved.
