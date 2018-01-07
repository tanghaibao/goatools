"""Used to operate on a sub-graph of a larger GO DAG."""

from collections import defaultdict

__copyright__ = "Copyright (C) 2016-2018, DV Klopfenstein, H Tang, All rights reserved."
__author__ = "DV Klopfenstein"

class GODagSmall(object):
    """Contains a sub-graph of the original obo from geneontology.org."""

    def __init__(self):
        # Sub-graph of input DAG. User has option of using a subset of children/parents
        self.go_sources = None
        self.go2obj = {}
        self.p_from_cs = defaultdict(set) # GO ids for child->parents
        self.c_from_ps = defaultdict(set) # GO ids for parent->children

    def num_goterms(self):
        """Return the number of GO terms in this GO DAG subgraph."""
        return len(self.go2obj)

    def get_edges(self):
        """Get the directed edges from GO term to GO term."""
        edge_from_to = []
        for parent, children in self.p_from_cs.items():
            for child in children:
                edge_from_to.append((child, parent))
        for parent, children in self.c_from_ps.items():
            for child in children:
                edge_from_to.append((child, parent))
        return edge_from_to
        
   
# Copyright (C) 2016-2018, DV Klopfenstein, H Tang, All rights reserved.

