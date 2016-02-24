"""Given GO ids and an obo, creates a small sub-graph DAG.

   Sub-graphs can be useful if we want to create shortcut paths and eliminate nodes.
"""

from collections import defaultdict
from goatools.godag_small import GODagSmall

__copyright__ = "Copyright (C) 2016, DV Klopfenstein, H Tang, All rights reserved."
__author__ = "DV Klopfenstein"

def get_godag(goids, obo_dag):
    """Given GO ids, return GODagSmall."""
    obj = OboToGoDagSmall(goids, obo_dag)
    return obj.godag

class OboToGoDagSmall(object):
    """Given a list of GO ids and the full obo dag, return a sub GO DAG graph."""

    def __init__(self, goids, obodag, **kws):
        self.obodag = obodag
        self.traverse_child = kws['traverse_child'] if 'traverse_child' in kws else False
        self.traverse_parent = kws['traverse_parent'] if 'traverse_parent' in kws else True
        self.seen_cids = set() if self.traverse_parent else None
        self.seen_pids = set() if self.traverse_child else None
        self.godag = GODagSmall()
        self._init(goids)
        assert self.traverse_child or self.traverse_parent, "NO EDGES IN GRAPH"

    def _init(self, goids):
        """Given GO ids and full obodag, create mini GO dag."""
        for goid in goids:
            goobj = self.obodag[goid]
            self.godag.go2obj[goid] = goobj
            # Traverse up parents
            if self.traverse_parent and goid not in self.seen_cids:
                self._traverse_parent_objs(goobj)
            # Traverse down children
            if self.traverse_child and goid not in self.seen_pids:
                self._traverse_child_objs(goobj)

    def _traverse_parent_objs(self, goobj_child):
        child_id = goobj_child.id
        # mark child as seen
        self.seen_cids.add(child_id)
        self.godag.go2obj[child_id] = goobj_child
        # Loop through parents of child object
        for parent_obj in goobj_child.parents:
            parent_id = parent_obj.id
            self.godag.p_from_cs[parent_id].add(child_id)
            # If parent has not been seen, traverse
            if parent_id not in self.seen_cids:
                self._traverse_parent_objs(parent_obj)
            
    def _traverse_child_objs(self, goobj_parent):
        parent_id = goobj_parent.id
        # mark parent as seen
        self.seen_pids.add(parent_id)
        self.godag.go2obj[parent_id] = goobj_parent
        # Loop through children
        for child_obj in goobj_parent.children:
            child_id = child_obj.id
            self.godag.p2cs[parent_id].add(child_id)
            # If child has not been seen
            if child_id not in self.seen_pids:
                self._traverse_child_objs(child_obj)

# Copyright (C) 2016, DV Klopfenstein, H Tang, All rights reserved.

