#!/usr/bin/env python3
"""Compare GOATOOLS Resnik scores and Yang Resnik scores"""

__copyright__ = "Copyright (C) 2019-2020, DV Klopfenstein, H Tang. All rights reserved."
__author__ = "DV Klopfenstein"

import os
from goatools.base import get_godag
from goatools.anno.gpad_reader import GpadReader
from goatools.associations import dnld_annotation
from goatools.semantic import TermCounts
from goatools.godag.consts import NS2NAMESPACE
from goatools.godag.consts import NS2GO


REPO = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..")


def test_termcnt_init():
    """Compare GOATOOLS Resnik scores and Yang Resnik scores"""
    godag = get_godag(os.path.join(REPO, 'go-basic.obo'))
    fin_gpad = os.path.join(REPO, 'goa_human.gpad')
    dnld_annotation(fin_gpad)

    # Load all annoations (BP, MF, CC)
    top_cnt_all = _run_full(fin_gpad, godag)

    # Load one annoation (BP, MF, CC) at a time
    top_cnt_ns = _run_each(fin_gpad, godag)

    # Compare different load methods
    assert top_cnt_all == top_cnt_ns


def _run_full(fin_gpad, godag):
    """Load all annoations (BP, MF, CC)"""
    annoobj = GpadReader(fin_gpad, godag=godag)
    id2gos = annoobj.get_id2gos('all')
    tcntobj = TermCounts(godag, id2gos)
    top_cnt_all = {}
    for nspc in ['BP', 'MF', 'CC']:
        top_ns = NS2GO[nspc]
        namespace = NS2NAMESPACE[nspc]
        top_cnt = tcntobj.gocnts[top_ns]
        top_cnt_all[nspc] = top_cnt
        assert top_cnt == tcntobj.aspect_counts[namespace]
    return top_cnt_all

def _run_each(fin_gpad, godag):
    """Load one annoation (BP, MF, CC) at a time"""
    top_cnt_ns = {}
    for nspc in ['BP', 'MF', 'CC']:
        namespace = NS2NAMESPACE[nspc]
        top_ns = NS2GO[nspc]
        annoobj = GpadReader(fin_gpad, godag=godag, namespaces={nspc})
        ns2assoc = annoobj.get_ns2assc()
        id2gos = ns2assoc[nspc]
        print('{NS} NUM ASSOC {N:6,}'.format(NS=nspc, N=len(ns2assoc[nspc])))
        tcntobj = TermCounts(godag, id2gos)
        print('{NS} GOATOOLS TERM COUNTS: {N}/{M}'.format(
            NS=nspc,
            N=tcntobj.gocnts[top_ns],
            M=tcntobj.aspect_counts[namespace]))
        top_cnt = tcntobj.gocnts[top_ns]
        top_cnt_ns[nspc] = top_cnt
        assert top_cnt == tcntobj.aspect_counts[namespace], '{NS} {A} != {B}'.format(
            NS=nspc, A=top_cnt, B=tcntobj.aspect_counts[namespace])
        assert top_cnt == max(tcntobj.gocnts.values())
    return top_cnt_ns


if __name__ == '__main__':
    test_termcnt_init()

# Copyright (C) 2019-2020, DV Klopfenstein, et al. All rights reserved.
