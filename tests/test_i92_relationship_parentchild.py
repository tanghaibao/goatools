#!/usr/bin/env python
"""Issue#92: parsing of part_of and other relationships as parent-children connections"""

from __future__ import print_function

import os
from goatools.base import get_godag

REPO = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..")

def test_i92():
    """Issue#92: parsing of part_of and other relationships as parent-children connections"""
    print('CWD', os.getcwd())
    # Loading heartjogging GO graph with the relationship tags
    fin_obo = "tests/data/heartjogging.obo"
    # Expected parents and upper(includes relationship)
    # GO:0007389: pattern specification process
    exp_par = {'GO:0007389': set(['GO:0008150', 'GO:0032501'])}
    exp_upper = {'GO:0007389':
                 set(['GO:0008150', 'GO:0032501', 'GO:0007275', 'GO:0048856', 'GO:0032502'])}
    # Expected children and upper(includes relationship_rev)
    # GO:0003143: embryonic heart tube morphogenesis
    exp_chi = {'GO:0003143': set(['GO:0003146'])}
    exp_lower = {'GO:0003143': set(['GO:0003146', 'GO:0003304'])}
    godag = get_godag(os.path.join(REPO, fin_obo), optional_attrs=['relationship'])
    # List GO terms that have relationships (relationship points to parent)
    _chk_parents(godag, exp_par, exp_upper)
    # List GO terms that have reverse relationships (relationship points to child)
    _chk_children(godag, exp_chi, exp_lower)


def _chk_children(godag, exp_chi, exp_lower):
    """Check all children(through 'is_a') and all lower(thru 'is_a' and all 'relationship_rev')."""
    print('\nrelationship_rev: Has relationship "children"')
    for goid, goterm in sorted(godag.items(), key=lambda t: t[1].depth):
        if goid == goterm.id:
            all_childids = sorted(goterm.get_all_children())
            all_lower = sorted(goterm.get_all_lower())
            # Check children through 'is_a'
            if goid in exp_chi:
                assert exp_chi[goid] == set(all_childids), "C[{}]: {}".format(goid, all_childids)
            # Check lower GO terms through 'is_a' and reverse 'relationship'
            if goid in exp_lower:
                assert exp_lower[goid] == set(all_lower)
            rel_revs = goterm.relationship_rev
            if rel_revs:
                print(goterm.namespace, goterm.depth, goterm.id, goterm.name)

def _chk_parents(godag, exp_par, exp_upper):
    """Check all parents(through 'is_a') and all upper(thru 'is_a' and all 'relationship')."""
    print('\nrelationship:')
    for goid, goterm in sorted(godag.items(), key=lambda t: [t[1].depth, t[1].id]):
        if goid == goterm.id:
            all_parentids = sorted(goterm.get_all_parents())
            all_upper = sorted(goterm.get_all_upper())
            if goterm.depth <= 3:
                print("\nget_all_parents[{GO}]: {NAME}".format(GO=goterm.id, NAME=goterm.name))
                print("get_all_parents[{GO}]: {GOs}".format(GO=goterm.id, GOs=all_parentids))
            # Check parents through 'is_a'
            if goid in exp_par:
                assert exp_par[goid] == set(all_parentids)
            # Check upper GO terms through 'is_a' and 'relationship'
            if goid in exp_upper:
                assert exp_upper[goid] == set(all_upper)
            # Print relationship
            rels = goterm.relationship
            if rels:
                print(goterm.namespace, goterm.depth, goterm.id, goterm.name)

# To create a plot to visually conform the test was set up correctly:
#
# $ go_plot -i test_gos.txt --obo=../goatools/tests/data/heartjogging.obo -r -o heartjogging_test.png

# FILE: test+gos.txt
# Contents
# #-------------------------------------------
# # Source for tests up GO hierachy
# GO:0007389#fffdda # Yellow
# 
# # Ancestors up hierarchy through is_a attribute
# GO:0008150#e0ffdc # Green
# GO:0032501#e0ffdc # Green
# 
# # Ancestors up hierarchy through relationship attributes
# GO:0007275#deddff # Purple 
# GO:0048856#deddff # Purple 
# GO:0032502#deddff # Purple
# 
# #-------------------------------------------
# # Source for tests down GO hierarchy
# GO:0003143#fffdda # Yellow  embryonic heart tube morphogenesis
# 
# # Ancestors up hierarchy through is_a attribute
# GO:0003146#e0ffdc # Green
# 
# # Ancestors up hierarchy through relationship attributes
# GO:0003304#deddff # Purple

if __name__ == '__main__':
    test_i92()
