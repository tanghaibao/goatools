"""Work with combinations of relationships, like *part_of* and *regulates*"""

from __future__ import print_function

__copyright__ = "Copyright (C) 2017-2019, DV Klopfenstein, H Tang et al. All rights reserved."

from goatools.godag.consts import RELATIONSHIP_SET


# pylint: disable=old-style-class
class RelationshipCombos:
    """Work with combinations of relationships, like *part_of* and *regulates*"""

    def __init__(self, godag):
        assert hasattr(next(iter(godag.values())), 'relationship'), "NO DAG RELATIONSHIPS"
        self.relationship_set = RELATIONSHIP_SET
        self.go2obj = {go:o for go, o in godag.items() if go == o.id}
        self.rels_all = self._init_rels_all()

    def get_relationship_combos(self):
        """Get all combinations of all lengths of relationship lists"""
        rels_combo = []
        num_rels = len(self.rels_all)
        print('GODAG relationships[{N}]: {Rs}'.format(N=num_rels, Rs=self.rels_all))
        for cnt in range(2**num_rels):
            idxs = [i for i, v in enumerate('{N:0{L}b}'.format(N=cnt, L=num_rels)) if v == '1']
            if idxs:
                rels_cur = set(self.rels_all[i] for i in idxs)
                rels_combo.append(rels_cur)
                # print('{N:0{L}b}'.format(N=cnt, L=num_rels), idxs, rels_cur)
        return rels_combo

    def chk_relationships_all(self):
        """Check that the list of relationships in consts is same as found in GODAG"""
        assert set(self.rels_all) == self.relationship_set, \
            set(self.rels_all).symmetric_difference(self.relationship_set)

    def _init_rels_all(self):
        """Get all relationships found in GODAG loaded with optional_attrs, relationship"""
        rels_all = set()
        for obj in self.go2obj.values():
            rels_all.update(obj.relationship.keys())
        return sorted(rels_all)


# Copyright (C) 2017-2019, DV Klopfenstein, H Tang et al. All rights reserved.
