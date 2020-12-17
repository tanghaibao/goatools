#!/usr/bin/env python3
"""Run code from issue #177, which is reporting a recursion error"""

from goatools.base import get_godag
from goatools.gosubdag.gosubdag import GoSubDag


def test_i177():
    """Run code from issue #177, which is reporting a recursion error"""
    go_id = 'GO:0050807'
    godag = get_godag('go.obo', optional_attrs='relationship')
    gosubdag_r0 = GoSubDag([go_id], godag, prt=None)
    print('{GO} ancestors: {P}'.format(
        GO=go_id,
        P=gosubdag_r0.rcntobj.go2ancestors[go_id]))


if __name__ == '__main__':
    test_i177()
