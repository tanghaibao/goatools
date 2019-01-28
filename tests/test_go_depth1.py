#!/usr/bin/env python
"""Test the use of CountRelatives to visualize the branche(s) of any GO term."""

import os
import sys
import collections as cx
from goatools.base import get_godag
from goatools.gosubdag.godag_rcnt import CountRelatives
#### from goatools.gosubdag.godag_depth1 import GoDepth1Letters
from goatools.gosubdag.rpt.wr_xlsx import GoDepth1LettersWr

REPO = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..")

def test_all(prt=sys.stdout):
    """Test initialization and operation of CountRelatives for GO term branch(s) visualization."""
    godag = get_godag(os.path.join(REPO, "go-basic.obo"), prt=sys.stdout)
    rcntobj = CountRelatives(godag)
    _wr_xlsx_d1(rcntobj)
    _run_get_letters_d1(rcntobj)
    _run_get_letters_d2(godag, rcntobj, prt)

def _run_get_letters_d1(rcntobj):
    """Test get_letter on goobjs at depth-01."""
    for goobj in rcntobj.depth2goobjs[1]:
        letters = rcntobj.get_parents_letters(goobj)
        assert len(letters) == 1
        assert letters[0] == rcntobj.goone2ntletter[goobj.id].D1

def _run_get_letters_d2(obo_dag, rcntobj, prt):
    """Print parent letters and descendant counts for all D1 and D2 GO Terms."""
    go2nt = __get_d1d2_goobjs(obo_dag, rcntobj)
    sortby = lambda t: [t[1].goobj.namespace, t[1].D1, -1*t[1].dcnt]
    for goid, ntlet in sorted(go2nt.items(), key=sortby):
        level = ntlet.goobj.level
        prt.write("{pre:1} {d:6,} {ABC:2} L{L:02} D{D:02} {GO} {NAME}\n".format(
            pre=_get_pre(level),
            ABC=ntlet.D1,
            d=ntlet.dcnt,
            L=level,
            D=ntlet.goobj.depth,
            GO=goid,
            NAME=ntlet.goobj.name))

def _get_pre(level):
    """Get prefix indicating a header GO or a user GO."""
    if level == 0:
        return ">"
    if level == 1:
        return "-"
    return ""

def __get_d1d2_goobjs(obo_dag, rcntobj):
    """Get goobjs at depth-00 (BP, MF, CC), depth-01 and depth-02."""
    go2nt = {}
    ntobj = cx.namedtuple("NtGoLetters", "D1 dcnt goobj")
    for goobj in sorted(obo_dag.values(), key=lambda nt: nt.depth):
        depth = goobj.depth
        if depth > 2:
            break
        if not goobj.is_obsolete:
            if goobj.id not in go2nt:
                go2nt[goobj.id] = ntobj(
                    D1=rcntobj.get_d1str(goobj),
                    dcnt=rcntobj.go2dcnt[goobj.id],
                    goobj=goobj)
    return go2nt

def _wr_xlsx_d1(rcntobj):
    """Test initialization."""
    d1wrobj = GoDepth1LettersWr(rcntobj)
    d1wrobj.prt_txt()
    d1wrobj.wr_xlsx("go_depth01_godag.xlsx")

if __name__ == '__main__':
    test_all()
