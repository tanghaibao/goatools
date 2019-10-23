#!/usr/bin/env python
"""Test that GoSubDag contains ancestors from only the user-specified relationships"""
# tests/test_gosubdag_relationships_i126.py
# goatools/gosubdag/gosubdag.py
# goatools/gosubdag/godag_rcnt.py
# goatools/gosubdag/godag_rcnt_init.py
# goatools/godag/go_tasks.py
# goatools/obo_parser.py

from __future__ import print_function

__copyright__ = "Copyright (C) 2016-2019, DV Klopfenstein, H Tang, All rights reserved."

import os
import sys
## import timeit
## import datetime
import collections as cx
from goatools.base import get_godag
from goatools.godag.consts import RELATIONSHIP_SET
from goatools.gosubdag.gosubdag import GoSubDag
from goatools.test_data.wr_subobo import WrSubObo


REPO = os.path.join(os.path.dirname(os.path.abspath(__file__)), "../")

# pylint: disable=line-too-long,unused-variable
def test_gosubdag_relationships(wr_new_obo_subset=False):
    """Test that GoSubDag contains ancestors from only the user-specified relationships"""

    # Leaf GO: viral triggering of virus induced gene silencing
    goid_chosen = 'GO:0060150'

    # Load GODag with all relationships
    fin_obo = os.path.join(REPO, "go-basic.obo")
    godag_r0 = get_godag(fin_obo, loading_bar=None)
    godag_r1 = get_godag(fin_obo, loading_bar=None, optional_attrs=['relationship'])

    file_sub = os.path.join(REPO, "tests/data/viral_gene_silence.obo")

    # Get all GO terms above this low-level GO ID using all relationships
    if wr_new_obo_subset:
        _wr_sub_obo(file_sub, goid_chosen, godag_r1, fin_obo)

    # RELATIONSHIPS: None
    gosubdag_r0 = GoSubDag(set([goid_chosen]), godag_r0)
    assert len(gosubdag_r0.rcntobj.go2ancestors[goid_chosen]) == 12

    # RELATIONSHIPS: ALL
    gosubdag_r1 = GoSubDag(set([goid_chosen]), godag_r1, relationships=True)
    assert gosubdag_r1.relationships == RELATIONSHIP_SET
        #### set(['part_of', 'regulates', 'positively_regulates', 'negatively_regulates'])
    assert len(gosubdag_r1.rcntobj.go2ancestors[goid_chosen]) == 50

    # RELATIONSHIPS: part_of
    gosubdag_rp = GoSubDag(set([goid_chosen]), godag_r1, relationships={'part_of'})
    assert gosubdag_rp.relationships == set(['part_of'])
    rp_par = gosubdag_rp.rcntobj.go2ancestors[goid_chosen]
    assert 'GO:0016441' not in gosubdag_rp.go2obj, '**FATAL: REGULATION TERM GoSubDag(part_of) go2obj'
    assert 'GO:0016441' not in rp_par, '**FATAL: REGULATION TERM GoSubDag(part_of) go2parents'

    # RELATIONSHIPS: regulates
    gosubdag_rr = GoSubDag(set([goid_chosen]), godag_r1, relationships={'regulates'})
    assert gosubdag_rr.relationships == set(['regulates'])
    rp_par = gosubdag_rr.rcntobj.go2ancestors[goid_chosen]
    # assert 'GO:0016441' not in gosubdag_rp.go2obj, '**FATAL: REGULATION TERM GoSubDag(part_of) go2obj'
    # assert 'GO:0016441' not in rp_par, '**FATAL: REGULATION TERM GoSubDag(part_of) go2parents'

    # RELATIONSHIPS: positively_regulates
    gosubdag_rp = GoSubDag(set([goid_chosen]), godag_r1, relationships={'positively_regulates'})
    assert gosubdag_rp.relationships == set(['positively_regulates'])
    rp_par = gosubdag_rp.rcntobj.go2ancestors[goid_chosen]

    # RELATIONSHIPS: negatively_regulates
    gosubdag_rn = GoSubDag(set([goid_chosen]), godag_r1, relationships={'negatively_regulates'})
    assert gosubdag_rn.relationships == set(['negatively_regulates'])
    rp_par = gosubdag_rn.rcntobj.go2ancestors[goid_chosen]

    # RELATIONSHIPS: regulates positively_regulates negatively_regulates
    regs = {'positively_regulates', 'negatively_regulates'}
    gosubdag_rnp = GoSubDag(set([goid_chosen]), godag_r1, relationships=regs)
    assert gosubdag_rnp.relationships == regs
    rp_par = gosubdag_rnp.rcntobj.go2ancestors[goid_chosen]

    _run_baseline_r0(gosubdag_r0, gosubdag_r1)

    # BASELINE r1: Test that GOTerm.get_all_upper() is the same as GoSubDag ancestors
    for goid, term in gosubdag_r1.go2obj.items():
        ancestors_r1 = gosubdag_r1.rcntobj.go2ancestors.get(goid, set())
        assert ancestors_r1 == term.get_all_upper()

    #### # Test that
    #### gosubdag_rp = GoSubDag(set([goid_chosen]), godag_r1, relationships={'part_of'}, prt=sys.stdout)
    #### for goid, dag_term in godag_r1.items():
    ####     if goid in gosubdag_r1.rcntobj.go2ancestors:
    ####         ancestors = gosubdag_rp.rcntobj.go2ancestors[goid]
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
    #### #assert 'GO:0016441' not in gosubdag_rp.rcntobj.go2ancestors['GO:0060150']
    #### assert 'GO:0016441' in gosubdag_r1.go2nt
    #### assert 'GO:0010467' in gosubdag_r1.go2nt


def _run_baseline_r0(gosubdag_r0, gosubdag_r1):
    """BASELINE r0: Test that GOTerm.get_all_parents() == GoSubDag ancestors"""
    r1_ancestors_more = set()
    # Loop through r0 GO IDs
    for goid, term in gosubdag_r0.go2obj.items():
        ancestors_r0 = gosubdag_r0.rcntobj.go2ancestors.get(goid, set())
        ancestors_r1 = gosubdag_r1.rcntobj.go2ancestors.get(goid, set())
        assert ancestors_r0 == term.get_all_parents()
        assert ancestors_r0.issubset(ancestors_r1)
        if len(ancestors_r0) < len(ancestors_r1):
            r1_ancestors_more.add(goid)
    assert len(r1_ancestors_more) != 0
    print('{N} r1 GO terms in GoSubDag have more ancestors than r0'.format(
        N=len(r1_ancestors_more)))
    # scripts/go_plot.py --go_file=i126_goids_baseline.txt -r --obo=tests/data/viral_gene_silence.obo -o i126_goids_baseline.png
    fout_gos = 'i126_goids_baseline.txt'
    with open(fout_gos, 'w') as prt:
        prt.write('#cafffb {SRC_GO}\n'.format(SRC_GO=next(iter(gosubdag_r0.go_sources))))
        _prt_goterms(r1_ancestors_more, gosubdag_r1.go2nt, prt)
        print('  WROTE: {GOs}'.format(GOs=fout_gos))

def _prt_goterms(goids, go2nt, prt):
    """Print details of GO terms"""
    fmt = ('#ffd1df {GO} # {NS} {dcnt:5} {childcnt:3} '
           'L{level:02} D{depth:02} R{reldepth:02} {D1:5} {REL} {rel} {GO_name}\n')
    nts = [nt for go, nt in go2nt.items() if go in goids]
    for ntd in sorted(nts, key=lambda nt: nt.dcnt, reverse=True):
        prt.write(fmt.format(**ntd._asdict()))

#cafffb GO:0060150
#ffd1df GO:0050794 # BP  8278  64 D03 R03 regulation of cellular process
#ffd1df GO:0019222 # BP  3382  20 D03 R03 regulation of metabolic process
#ffd1df GO:0048522 # BP  2417  65 D04 R04 positive regulation of cellular process
#ffd1df GO:0060255 # BP  2130  20 D04 R04 regulation of macromolecule metabolic process
#ffd1df GO:0010468 # BP   862  20 D05 R05 regulation of gene expression
#ffd1df GO:0060968 # BP    53   4 D06 R08 regulation of gene silencing
#ffd1df GO:0060147 # BP    24   4 D07 R09 regulation of posttranscriptional gene silencing
#ffd1df GO:0060148 # BP     8   3 D08 R10 positive regulation of posttranscriptional gene silencing
#ffd1df GO:0060150 # BP     0   0 D09 R11  viral triggering of virus induced gene silencing

# - Generate GO DAG subset for this test ---------------------------------------------------------
def _wr_sub_obo(fout_obo, goid_chosen, godag_r1, fin_obo):
    """Sub plot used for visualizing this test file's elements"""
    # Load GO-DAG: Load optional 'relationship'
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
    # Plot obo subset
    pat_r1 = '{REPO}/scripts/go_plot.py {GO} -o {PNG} -r'
    pat_r0 = '{REPO}/scripts/go_plot.py {GO} -o {PNG}'
    os.system(pat_r1.format(REPO=REPO, PNG=fout_obo.replace('.obo', '_r1.png'), GO=goid_chosen))
    os.system(pat_r0.format(REPO=REPO, PNG=fout_obo.replace('.obo', '_r0.png'), GO=goid_chosen))

def _get_leafs_w_relsinhier(rels_usr, gosubdag_r1):
    """Get GO IDs that have all relationships up their hierarchy."""
    gos_r1_relsinhier = set()
    goids_leaf = set(o.id for o in gosubdag_r1.go2obj.values() if not o.children)
    for goid in goids_leaf:
        go_parents = gosubdag_r1.rcntobj.go2ancestors[goid]
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
