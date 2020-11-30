"""Wang's termwise semantic similarity for GO terms"""
# http://127.0.0.1:31762/library/GOSemSim/doc/GOSemSim.html
# https://cran.r-project.org/doc/manuals/r-release/R-intro.html#An-introductory-session
# TOP 20179076 R. HA..c  98 4 2010   332  1  13 au[06](Guangchuang Yu)
#     GOSemSim: an R package for measuring semantic similarity among GO terms and gene products
# https://github.com/mojaie/pygosemsim
# https://bioconductor.org/packages/release/bioc/html/GOSemSim.html

__copyright__ = "Copyright (C) 2020-present, DV Klopfenstein. All rights reserved."
__author__ = "DV Klopfenstein"

##import timeit
##from goatools.godag.prttime import prt_hms

from sys import stdout
from goatools.semsim.termwise.dag_a import DagA
from goatools.godag.go_tasks import get_go2ancestors
from goatools.godag.go_tasks import get_go2depth

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
        self.rels = self._init_rels(relationships)
        self.w_e = self._init_edge_weight_factor(rel2scf)
        # Get GO Wang DAGs
        self.go2dag = self._init_go2dag(goids)

    def get_sim(self, go_a, go_b):
        """Get Wang's semantic similarity between two GO terms"""
        if self._not_loaded(go_a, go_b):
            return None
        dag_a = self.go2dag[go_a]
        dag_b = self.go2dag[go_b]
        goids_ab = set(dag_a.goids).intersection(dag_b.goids)
        s_a = dag_a.get_svalues(goids_ab)
        s_b = dag_b.get_svalues(goids_ab)
        s_ab = sum([*s_a, *s_b])
        return s_ab/(dag_a.get_sv() + dag_b.get_sv())

    def prt_cfg(self, prt=stdout):
        """Print reseacher-specified Wang configuration"""
        prt.write('Wang Semantic Similarity Configuration:\n')
        prt.write('    Optional relationships: {Rs}\n'.format(Rs=' '.join(sorted(self.rels))))
        prt.write('    Edge weights:\n')
        for rel, weight in self.w_e.items():
            prt.write('        {W:.8f} {REL}\n'.format(W=weight, REL=rel))
        prt.write('\n')

    def _not_loaded(self, go_a, go_b):
        """Check that both GO IDs are in the go2subdag dict"""
        if go_a not in self.go2dag:
            print('**ERROR: {GO} NOT LOADED INTO SsWang'.format(GO=go_a))
            return True
        if go_b not in self.go2dag:
            print('**ERROR: {GO} NOT LOADED INTO SsWang'.format(GO=go_b))
            return True
        return False

    def _init_edge_weight_factor(self, rel2scf):
        """Initialize semantic contribution factor (scf) for weights for edge types (w_e)"""
        d_rel2scf = self.dflt_rel2scf
        rels = self.rels
        # If there are no user-provided edge weights
        ret = {'is_a': d_rel2scf['is_a']}
        if not rel2scf:
            if not rels:
                return ret
            for rel in rels:
                ret[rel] = d_rel2scf[rel]
            return ret
        # If there are user-provided edge weights
        if 'is_a' in rel2scf:
            ret = {'is_a': rel2scf['is_a']}
        if not rels:
            return ret
        for rel in rels:
            ret[rel] = rel2scf[rel] if rel in rel2scf else d_rel2scf[rel]
        return ret

    def _init_go2dag(self, goids):
        """Get all GO IDs in the DAG above and including GO IDs in goids arg"""
        # GO terms provided by user
        ##tic = timeit.default_timer()
        s_godag = self.godag
        rels = self.rels
        go_set_all = set(goids)
        go_set_cur = go_set_all.intersection(s_godag.keys())
        if go_set_cur != go_set_all:
            self._go_not_found(go_set_cur, go_set_all)
        # Ancestor GO terms for each user GO term
        ##tic = prt_hms(tic, '_init_go2dag GO IDs not found')
        go2ancestors = get_go2ancestors(self._get_goobjs(go_set_cur), rels)
        ##tic = prt_hms(tic, '_init_go2dag go2ancestors')
        go2depth = self._get_go2depth(go2ancestors, rels)
        ##tic = prt_hms(tic, '_init_go2dag go2depth')
        w_e = self.w_e
        # pylint: disable=line-too-long
        go2dag = {go:DagA(go, ancestors, go2depth, w_e, s_godag) for go, ancestors in go2ancestors.items()}
        ##tic = prt_hms(tic, '_init_go2dag DagA')
        # Add alt GO IDs
        for go_alt in go_set_cur.difference(go2ancestors.keys()):
            go_term = s_godag[go_alt]
            go_main = go_term.item_id
            go_depth = go_term.depth
            if go_depth != 0:
                go2dag[go_alt] = go2dag[go_main]
            elif go_depth == 0:
                go2dag[go_alt] = DagA(go_main, {}, go2depth, w_e, s_godag)
        ##tic = prt_hms(tic, '_init_go2dag ALT GO IDs')
        return go2dag

    def _get_goobjs(self, goids):
        """Get ancestors of GO IDs"""
        s_godag = self.godag
        return {s_godag[go] for go in goids}

    def _get_go2depth(self, go2ancestors, rels):
        """Get depth of all GO IDs and their ancestors"""
        goids_all = set.union(*go2ancestors.values(), set(go2ancestors.keys()))
        goobjs = self._get_goobjs(goids_all)
        return get_go2depth(goobjs, rels)

    @staticmethod
    def _go_not_found(go_set_cur, go_set_all):
        """GO IDs provided by researcher which are not found"""
        print('**WARNING: GO IDs NOT FOUND: {GOs}'.format(
            GOs=sorted(go_set_all.difference(go_set_cur))))

    def _init_rels(self, relationships):
        """Return valid relationships are valid"""
        if relationships is None:
            return {}
        if not hasattr(next(iter(self.godag.values())), 'relationship'):
            raise RuntimeError('**ERROR: SsWang GODag not loaded with relationships')
        rels = set(relationships)
        if rels.issubset(self.dflt_rel2scf.keys()):
            return rels
        print('**ERROR: INVALID RELATIONSHIPS: {Rs}'.format(
            Rs=rels.difference(self.dflt_rel2scf.keys())))
        return rels.intersection(self.dflt_rel2scf.keys())


# Copyright (C) 2020-present DV Klopfenstein. All rights reserved.
