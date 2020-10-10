"""Work with combinations of relationships, like *part_of* and *regulates*"""

from __future__ import print_function

__copyright__ = "Copyright (C) 2017-present, DV Klopfenstein, H Tang et al. All rights reserved."

from goatools.godag.consts import RELATIONSHIP_SET


class RelationshipCombos:
    """Work with combinations of relationships, like *part_of* and *regulates*"""

    def __init__(self, godag):
        self.godag_rels_loaded = hasattr(next(iter(godag.values())), 'relationship')
        self.relationship_set = RELATIONSHIP_SET
        self.go2obj = {go:o for go, o in godag.items() if go == o.id}
        self.dag_rels = self._init_dag_relationships()

    def get_relationship_combos(self):
        """Get all combinations of all lengths of relationship lists"""
        rels_combo = []
        num_rels = len(self.dag_rels)
        dag_rels = sorted(self.dag_rels)
        print('GODAG relationships[{N}]: {Rs}'.format(N=num_rels, Rs=dag_rels))
        for cnt in range(2**num_rels):
            idxs = [i for i, v in enumerate('{N:0{L}b}'.format(N=cnt, L=num_rels)) if v == '1']
            if idxs:
                rels_cur = set(dag_rels[i] for i in idxs)
                rels_combo.append(rels_cur)
                # print('{N:0{L}b}'.format(N=cnt, L=num_rels), idxs, rels_cur)
        return rels_combo

    def chk_relationships_all(self):
        """Check that the list of relationships in consts is same as found in GODAG"""
        assert set(self.dag_rels) == self.relationship_set, \
            set(self.dag_rels).symmetric_difference(self.relationship_set)

    def get_set(self, relationships_arg):
        """Return a set of relationships found in all subset GO Terms."""
        if relationships_arg:
            if self.godag_rels_loaded:
                relationships_dag = self.dag_rels
                if relationships_arg is True:
                    return relationships_dag
                relationshipset_usr = self._get_set_rel(relationships_arg)
                if relationshipset_usr:
                    self._chk_expected(relationshipset_usr, relationships_dag)
                    return relationships_dag.intersection(relationshipset_usr)
                print('**WARNING: UNKNOWN GODag relationships({R}). EXPECTED VALUES: {Rs}'.format(
                    R=relationships_arg, Rs=' '.join(sorted(relationships_dag))))
            else:
                err = ("""**WARNING: IGNORING(relationships={R}); NO GODag RELATIONSHIPS LOADED """
                       """W/ARG, optional_attrs=['relationship']""")
                print(err.format(R=relationships_arg))
        return set()

    @staticmethod
    def _get_set_rel(relationships):
        """Return a set containing with prospective relationships"""
        if isinstance(relationships, set):
            return relationships
        if isinstance(relationships, list):
            return set(relationships)
        if isinstance(relationships, str):
            return set([relationships,])
        return None

    @staticmethod
    def _chk_expected(relationships_usr, relationships_dag):
        """Check that user relationships were found"""
        rels_unexpected = relationships_usr.difference(relationships_dag)
        if not rels_unexpected:
            return
        print('**NOTE: RELATIONSHIPS IN GODag: SEEN({Rs}) NOT_SEEN({R})'.format(
            R=' '.join(sorted(rels_unexpected)), Rs=' '.join(sorted(relationships_dag))))

    def _init_dag_relationships(self):
        """Return all relationships seen in GO Dag subset."""
        relationship_set = set()
        if not self.godag_rels_loaded:
            return relationship_set
        for goterm in self.go2obj.values():
            if goterm.relationship:
                relationship_set.update(goterm.relationship)
            if goterm.relationship_rev:
                relationship_set.update(goterm.relationship_rev)
        return relationship_set

# Copyright (C) 2017-present, DV Klopfenstein, H Tang et al. All rights reserved.
