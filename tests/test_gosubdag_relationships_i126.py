#!/usr/bin/env python
"""Plot both the standard 'is_a' field and the optional 'part_of' relationship."""

from __future__ import print_function

__copyright__ = "Copyright (C) 2016-2019, DV Klopfenstein, H Tang, All rights reserved."

import os
import sys
## import timeit
## import datetime
import collections as cx
from goatools.base import get_godag
from goatools.gosubdag.gosubdag import GoSubDag
from goatools.test_data.wr_subobo import WrSubObo


REPO = os.path.join(os.path.dirname(os.path.abspath(__file__)), "../")

# pylint: disable=line-too-long,unused-variable
def test_gosubdag_relationships(wr_new_obo_subset=False):
    """Plot both the standard 'is_a' field and the 'part_of' relationship."""

    goid_chosen = 'GO:0060150'

    # Load GODag with all relationships
    fin_obo = os.path.join(REPO, "go-basic.obo")
    # godag_r0 = get_godag(file_obo, loading_bar=None)
    godag_r1 = get_godag(fin_obo, loading_bar=None, optional_attrs=['relationship'])

    file_sub = os.path.join(REPO, "tests/data/viral_gene_silence.obo")

    # Get all GO terms above this low-level GO ID using all relationships
    if wr_new_obo_subset:
        _wr_sub_obo(file_sub, goid_chosen, godag_r1, fin_obo)

    gosubdag_r0 = GoSubDag(set([goid_chosen]), godag_r1, relationships=False, prt=sys.stdout)
    gosubdag_r1 = GoSubDag(set([goid_chosen]), godag_r1, relationships=True, prt=sys.stdout)

    _run_baseline_r0(gosubdag_r0, gosubdag_r1)

    # BASELINE r1: Test that GOTerm.get_all_upper() is the same as GoSubDag ancestors
    for goid, term in gosubdag_r1.go2obj.items():
        ancestors_r1 = gosubdag_r1.rcntobj.go2parents[goid]
        assert ancestors_r1 == term.get_all_upper()

    #### # Test that
    #### gosubdag_rp = GoSubDag(set([goid_chosen]), godag_r1, relationships={'part_of'}, prt=sys.stdout)
    #### for goid, dag_term in godag_r1.items():
    ####     if goid in gosubdag_r1.rcntobj.go2parents:
    ####         ancestors = gosubdag_rp.rcntobj.go2parents[goid]
    ####         sub_term = gosubdag_rp.go2obj[goid]
    ####         reldict = sub_term.relationship.items()
    ####         # print(goid)
    ####         # print('DAG', sorted(dag_term.get_all_upper()))
    ####         # print('SUB', sorted(sub_term.get_all_upper()))
    ####         # print('ANS', sorted(ancestors))
    ####         # for rel, pterms in cx.OrderedDict(reldict).items():
    ####         #     print(rel, ' '.join(sorted(o.id for o in pterms)))
    ####         # print('')
    #### print(gosubdag_rp.relationships)
    #### #assert 'GO:0016441' not in gosubdag_rp.rcntobj.go2parents['GO:0060150']
    #### assert 'GO:0016441' in gosubdag_r1.go2nt
    #### assert 'GO:0010467' in gosubdag_r1.go2nt


def _run_baseline_r0(gosubdag_r0, gosubdag_r1):
    # BASELINE r0: Test that GOTerm.get_all_parents() is the same as GoSubDag ancestors
    r1_ancestors_more = set()
    for goid, term in gosubdag_r0.go2obj.items():
        ancestors_r0 = gosubdag_r0.rcntobj.go2parents[goid]
        ancestors_r1 = gosubdag_r1.rcntobj.go2parents[goid]
        assert ancestors_r0 == term.get_all_parents()
        assert ancestors_r0.issubset(ancestors_r1)
        if len(ancestors_r0) < len(ancestors_r1):
            r1_ancestors_more.add(goid)
    assert len(r1_ancestors_more) != 0
    _prt_goterms(r1_ancestors_more, gosubdag_r1.go2nt)

def _prt_goterms(goids, go2nt):
    """Print details of GO terms"""
    print('{N} r1 GO terms in GoSubDag have more ancestors than r0'.format(N=len(goids)))
    fmt = '{GO} # {NS} {dcnt:5} {childcnt:3} L{level:02} D{depth:02} R{reldepth:02} {D1:5} {REL} {rel} {GO_name}'
    nts = [nt for go, nt in go2nt.items() if go in goids]
    for ntd in sorted(nts, key=lambda nt: nt.dcnt, reverse=True):
        print(fmt.format(**ntd._asdict()))


# - Generate GO DAG subset for this test ---------------------------------------------------------
def _wr_sub_obo(fout_obo, goid_chosen, godag_r1, fin_obo):
    """Sub plot used for visualizing this test file's elements"""
    # Load GO-DAG: Load optional 'relationship'
    pngpat = '{REPO}/scripts/go_plot.py -r {GO} -o {PNG}'
    godag = {go:o for go, o in godag_r1.items() if go == o.item_id}
    _prt_rtel_ctr(godag)
    rels_all = set(['part_of', 'regulates', 'negatively_regulates', 'positively_regulates'])
    goids_leaf_all = set(o.id for o in godag.values() if not o.children)
    gosubdag_r1 = GoSubDag(goids_leaf_all, godag, relationships=True, prt=sys.stdout)
    goids_src_r1_all = _get_leafs_w_relsinhier(rels_all, gosubdag_r1)
    gosubdag_r1.prt_goids(goids_src_r1_all)
    # Pick one of the GO IDs as a source for the subset DAG
    gosubdag_viral = GoSubDag({goid_chosen}, godag, relationships=True, prt=sys.stdout)
    goids_viral = set(gosubdag_viral.go2obj.keys())
    with open(fout_obo, 'w') as prt:
        WrSubObo.prt_goterms(fin_obo, goids_viral, prt)
        print('{N} GO IDs WROTE: {OBO}'.format(N=len(goids_viral), OBO=fout_obo))
    os.system(pngpat.format(REPO=REPO, PNG=fout_obo.replace('obo', 'png'), GO=goid_chosen))

def _get_leafs_w_relsinhier(rels_usr, gosubdag_r1):
    """Get GO IDs that have all relationships up their hierarchy."""
    gos_r1_relsinhier = set()
    goids_leaf = set(o.id for o in gosubdag_r1.go2obj.values() if not o.children)
    for goid in goids_leaf:
        go_parents = gosubdag_r1.rcntobj.go2parents[goid]
        rels = set(k for p in go_parents for k in gosubdag_r1.go2obj[p].relationship.keys())
        if rels == rels_usr:
            gos_r1_relsinhier.add(goid)
    return gos_r1_relsinhier

def _prt_rtel_ctr(godag):
    """Print the count of relationships."""
    objs_r1_all = set(o for o in godag.values() if o.relationship.keys())
    octr = cx.Counter(k for o in objs_r1_all for k in o.relationship.keys())
    # objs_r1_sub = set(o.id for o in objs_r1_all if not rels_all.isdisjoint(o.relationship.keys()))
    print('{N:6,} GO Terms have relationships.'.format(N=len(objs_r1_all)))
    for key, cnt in octr.most_common():
        print('{N:6,} {REL}'.format(N=cnt, REL=key))

# def _chk_child_parent(go2o_dag, go2o_sub):
#     """Check the differences between the two go2obb dicts."""
#     pass

if __name__ == '__main__':
    test_gosubdag_relationships(len(sys.argv) != 1)

# Copyright (C) 2016-2019, DV Klopfenstein, H Tang, All rights reserved.
