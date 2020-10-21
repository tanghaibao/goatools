"""Wang's termwise semantic similarity for GO terms"""

__copyright__ = "Copyright (C) 2020-present, DV Klopfenstein. All rights reserved."
__author__ = "DV Klopfenstein"

from goatools.gosubdag.gosubdag import GoSubDag
from goatools.semsim.termwise.dag_a import DagA
from goatools.godag.go_tasks import get_go2ancestors


class SsWang:
    """Wang's termwise semantic similarity for GO terms"""

    # Default semantic contribution factor (scf) for weights for edge types (w_e)
    dflt_rel2scf = {
        'is_a': 0.8,
        'part_of': 0.6,
        'regulates': 0.6,
        'negatively_regulates': 0.6,
        'positively_regulates': 0.6,
    }

    def __init__(self, goids, godag, relationships=None, rel2scf=None):
        self.godag = godag
        self.rels = relationships
        self.w_e = self._init_edge_weight_factor(rel2scf)
        self.go2subdag = {}
        self.go2ancestors = self._init_go2ancestors(goids, godag, relationships)

    def add_goid(self, goid, prt=None):
        """Add a GO ID which will be compared using semantic similarity"""
        self.add_goids([goid], prt)

    def add_goids(self, goids, prt=None):
        """Add GO IDs which will be compared using semantic similarity"""
        # go2svalue = wang.get_go2svalue('go:0043231')
        s_godag = self.godag
        s_rels = self.rels
        s_go2subdag = self.go2subdag
        s_rel2scf = self.w_e
        for goid in goids:
            if goid in s_godag:
                gosubdag = GoSubDag([goid], s_godag, s_rels, prt=prt)
                dag = DagA(goid, gosubdag, s_rel2scf)
                s_go2subdag[goid] = dag
            else:
                print('**WARNING: {GO} NOT FOUND'.format(GO=goid))

    def get_sim(self, go_a, go_b):
        """Get Wang's semantic similarity between two GO terms"""
        if self._not_loaded(go_a, go_b):
            return None
        dag_a = self.go2subdag[go_a]
        dag_b = self.go2subdag[go_b]
        gos_ab = set(dag_a.goids).intersection(dag_b.goids)
        s_a = dag_a.get_svalues(gos_ab)
        s_b = dag_b.get_svalues(gos_ab)
        s_ab = sum([*s_a, *s_b])
        return s_ab/(dag_a.get_sv() + dag_b.get_sv())

    def _not_loaded(self, go_a, go_b):
        """Check that both GO IDs are in the go2subdag dict"""
        if go_a not in self.go2subdag:
            print('**ERROR: {GO} NOT LOADED INTO SsWang'.format(GO=go_a))
            return True
        if go_b not in self.go2subdag:
            print('**ERROR: {GO} NOT LOADED INTO SsWang'.format(GO=go_b))
            return True
        return False

    def _init_edge_weight_factor(self, rel2scf):
        """Initialize semantic contribution factor (scf) for weights for edge types (w_e)"""
        if rel2scf is None:
            return self.dflt_rel2scf
        ret = dict(self.dflt_rel2scf)
        rels_exp = set(self.dflt_rel2scf.keys())
        for rel, val in rel2scf.items():
            if rel in rels_exp:
                ret[rel] = val
            else:
                print('**ERROR: UNEXPECTED SEMANTIC CONTRIBUTION FACTOR: {K} = {V}'.format(
                    K=rel, V=val))
        return ret

    def _init_go2ancestors(self, goids, godag, relationships):
        go_set_all = set(goids)
        go_set_cur = go_set_all.intersection(godag.keys())
        if go_set_cur != go_set_all:
            self._go_not_found(go_set_cur, go_set_all)
        goobjs = [godag[go] for go in go_set_cur]
        return get_go2ancestors(goobjs, relationships)

    @staticmethod
    def _go_not_found(go_set_cur, go_set_all):
        """GO IDs provided by researcher which are not found"""
        print('**WARNING: GO IDs NOT FOUND: {GOs}'.format(
            GOs=sorted(go_set_all.difference(go_set_cur))))


# Copyright (C) 2020-present DV Klopfenstein. All rights reserved.
