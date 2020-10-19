"""A GO term, A, can be represented as DAG_a = (A, T_a, E_a), aka a GoSubDag"""

__copyright__ = "Copyright (C) 2020-present, DV Klopfenstein. All rights reserved."
__author__ = "DV Klopfenstein"


class DagA:
    """A GO term, A, can be represented as DAG_a = (A, T_a, E_a), aka a GoSubDag"""

    def __init__(self, go_a, gosubdag, rel2scf):
        self.go_usr = go_a
        self.go_a = gosubdag.go2obj[go_a].item_id
        self.gosubdag = gosubdag
        self.go2svalue = self._init_go2svalue(rel2scf)

    def get_sv(self):
        """Get the semantic value of GO term A"""
        return sum(self.go2svalue.values())

    def get_svalues(self, goids):
        """Get svalues for given IDs"""
        s_go2svalue = self.go2svalue
        return [s_go2svalue[go] for go in goids]

    def _init_go2svalue(self, rel2scf):
        """S-value: the contribution of GO term, t, to the semantics of GO term, A"""
        go2svalue = {self.go_a: 1.0}
        s_go2obj = self.gosubdag.go2obj
        s_rels = self.gosubdag.relationships
        terms_a = set(self.gosubdag.rcntobj.go2ancestors[self.go_a])
        ancestors_sorted = self._get_sorted(terms_a)
        terms_a.add(self.go_a)
        for ancestor_id, _ in ancestors_sorted:
            ## print(ancestor_id, ntd)
            goterm = s_go2obj[ancestor_id]
            svals = []
            weight = rel2scf['is_a']
            for cobj in goterm.children:
                if cobj.item_id in terms_a:
                    svals.append(weight*go2svalue[cobj.item_id])
            weight = rel2scf['part_of']
            for rel in s_rels:
                if rel in goterm.relationship_rev:
                    for cobj in goterm.relationship_rev[rel]:
                        if cobj.item_id in terms_a:
                            svals.append(weight*go2svalue[cobj.item_id])
            if svals:
                go2svalue[ancestor_id] = max(svals)
            ## print(ancestor_id, max(svals))
        return go2svalue

    def _get_sorted(self, ancestors):
        """Get the sorted ancestors"""
        go2nt = self.gosubdag.get_go2nt(ancestors)
        if self.gosubdag.relationships:
            return sorted(go2nt.items(), key=lambda t: t[1].reldepth, reverse=True)
        return sorted(go2nt.items(), key=lambda t: t[1].depth, reverse=True)


# Copyright (C) 2020-present, DV Klopfenstein. All rights reserved."
