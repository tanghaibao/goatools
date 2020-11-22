"""A GO term, A, can be represented as DAG_a = (A, T_a, E_a), aka a GoSubDag"""

__copyright__ = "Copyright (C) 2020-present, DV Klopfenstein. All rights reserved."
__author__ = "DV Klopfenstein"

## import timeit
## from goatools.godag.prttime import prt_hms


class DagA:
    """A GO term, A, can be represented as DAG_a = (A, T_a, E_a), aka a GoSubDag"""

    def __init__(self, go_a, ancestors, go2depth, w_e, godag):
        self.go_a = go_a
        self.ancestors = ancestors
        #tic = timeit.default_timer()
        self.goids = self._init_goids()
        #prt_hms(tic, '\nDagA INIT GO IDs')
        self.go2svalue = self._init_go2svalue(go2depth, w_e, godag)
        #prt_hms(tic, 'DagA SVALUES')

    def get_sv(self):
        """Get the semantic value of GO term A"""
        return sum(self.go2svalue.values())

    def get_svalues(self, goids):
        """Get svalues for given IDs"""
        s_go2svalue = self.go2svalue
        return [s_go2svalue[go] for go in goids]

    def _init_go2svalue(self, go2depth, w_e, godag):
        """S-value: the contribution of GO term, t, to the semantics of GO term, A"""
        #tic = timeit.default_timer()
        go2svalue = {self.go_a: 1.0}
        if not self.ancestors:
            return go2svalue
        terms_a = self.goids
        w_r = {r:v for r, v in w_e.items() if r != 'is_a'}
        #prt_hms(tic, 'DagA edge weights wo/is_a')
        for ancestor_id in self._get_sorted(go2depth):
            goterm = godag[ancestor_id]
            weight = w_e['is_a']
            svals = [weight*go2svalue[o.item_id] for o in goterm.children if o.item_id in terms_a]
            for rel, weight in w_r.items():
                if rel in goterm.relationship_rev:
                    for cobj in goterm.relationship_rev[rel]:
                        if cobj.item_id in terms_a:
                            svals.append(weight*go2svalue[cobj.item_id])
            if svals:
                go2svalue[ancestor_id] = max(svals)
            ## print(ancestor_id, max(svals))
        return go2svalue

    def _get_sorted(self, go2depth):
        """Get the sorted ancestors"""
        #tic = timeit.default_timer()
        go2dep = {go:go2depth[go] for go in self.ancestors}
        go_dep = sorted(go2dep.items(), key=lambda t: t[1], reverse=True)
        gos, _ = zip(*go_dep)
        #prt_hms(tic, 'DagA SORTED')
        return gos

    def _init_goids(self):
        """Return all GO IDs in GO_a's GODAG"""
        goids = set(self.ancestors)
        goids.add(self.go_a)
        return goids


# Copyright (C) 2020-present, DV Klopfenstein. All rights reserved."
