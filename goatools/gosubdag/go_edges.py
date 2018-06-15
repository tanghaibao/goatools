"""Creates and manages edges from one GO term to another GO term."""

__copyright__ = "Copyright (C) 2016-2018, DV Klopfenstein, H Tang, All rights reserved."
__author__ = "DV Klopfenstein"

from collections import defaultdict


def get_edgesobj(gosubdag, **kws):
    """Return specfied GoSubDag initialization object."""
    # Keyword args (kws):
    #     1. dst_srcs_list  Used for edges pruned such that only GO terms
    #                       are retained which are between the sets of dst & srcs.
    #     2  traverse_parent & traverse_child
    #                       Used to generate a GoSubDag with all parent terms and/or
    #                       all child terms, without pruning any paths.
    # Call function, get_edgesobj, with:
    #     get_edgesobj(go2obj, dst_srcs_list=...)
    # Or any of:
    #     get_edgesobj(go2obj, go_sources=...)
    #     get_edgesobj(go2obj, go_sources=..., traverse_parent=...,)
    #     get_edgesobj(go2obj, go_sources=..., traverse_child=...,)
    #     get_edgesobj(go2obj, go_sources=..., traverse_parent=..., traverse_child=...,)
    edgeobj = _get_edgesobj(gosubdag, **kws)
    rm_gos = kws.get('rm_gos')
    if rm_gos is not None:
        edgeobj.rm_gos(rm_gos)
    return edgeobj

def _get_edgesobj(gosubdag, **kws):
    """Return specfied GoSubDag initialization object."""
    # Keyword args (kws):
    #     1. dst_srcs_list  Used for edges pruned such that only GO terms
    #                       are retained which are between the sets of dst & srcs.
    #     2  traverse_parent & traverse_child
    #                       Used to generate a GoSubDag with all parent terms and/or
    #                       all child terms, without pruning any paths.
    # Call function, get_edgesobj, with:
    #     get_edgesobj(go2obj, dst_srcs_list=...)
    # Or any of:
    #     get_edgesobj(go2obj, go_sources=...)
    #     get_edgesobj(go2obj, go_sources=..., traverse_parent=...,)
    #     get_edgesobj(go2obj, go_sources=..., traverse_child=...,)
    #     get_edgesobj(go2obj, go_sources=..., traverse_parent=..., traverse_child=...,)
    dst_srcs_list = kws.get('dst_srcs_list', None)
    if dst_srcs_list is not None:
        return EdgesPath(gosubdag, dst_srcs_list)
    return EdgesRelatives(gosubdag,
                          kws.get('traverse_parent', True),
                          kws.get('traverse_child', False))

# -- Base Class ----------------------------------------------------------------
class EdgesBase(object):
    """Base class for GoEdges class."""

    def __init__(self, gosubdag):
        self.gosubdag = gosubdag
        self.go2obj = gosubdag.go2obj
        self.relationships = gosubdag.relationships
        # Set by derived edge class
        self.edges = []  # Lists of (goid_child, goid_parent)
        self.edges_rel = {}

    def rm_gos(self, rm_goids):
        """Remove any edges that contain user-specified edges."""
        self.edges = self._rm_gos_edges(rm_goids, self.edges)
        self.edges_rel = self._rm_gos_edges_rel(rm_goids, self.edges_rel)

    def _rm_gos_edges_rel(self, rm_goids, edges_rel):
        """Remove any relationship that contain user-specified edges."""
        edges_ret = {}
        for rname, edges_cur in edges_rel.items():
            edges_new = self._rm_gos_edges(rm_goids, edges_cur)
            if edges_new:
                edges_ret[rname] = edges_new
        return edges_ret

    @staticmethod
    def _rm_gos_edges(rm_goids, edges_all):
        """Remove any is_a edges that contain user-specified edges."""
        edges_reduced = []
        for goid_child, goid_parent in sorted(edges_all, key=lambda t: t[1]):
            if goid_child not in rm_goids and goid_parent not in rm_goids:
                edges_reduced.append((goid_child, goid_parent))
        return edges_reduced

    def get_all_edge_nodes(self):
        """Return a list of all GO IDs that are connected to edges."""
        edge_nodes = set(e for es in self.edges for e in es)
        for edges in self.edges_rel.values():
            rel_nodes = set(e for es in edges for e in es)
            edge_nodes.update(rel_nodes)
        return edge_nodes

    def chk_edges(self):
        """Check that all edge nodes exist in local subset."""
        goids = set(self.go2obj)
        self.chk_edges_nodes(self.edges, goids, "is_a")
        for reltype, edges in self.edges_rel.items():
            self.chk_edges_nodes(edges, goids, reltype)

    @staticmethod
    def chk_edges_nodes(edges, nodes, name):
        """Check that user specified edges have a node which exists."""
        edge_nodes = set(e for es in edges for e in es)
        missing_nodes = edge_nodes.difference(nodes)
        assert not missing_nodes, "MISSING: {GOs}\n{NM} EDGES MISSING {N} NODES (OF {T})".format(
            NM=name, N=len(missing_nodes), T=len(edge_nodes), GOs=missing_nodes)

    def get_c2ps(self):
        """Set child2parents dict for all parents used in this set of edges."""
        c2ps = defaultdict(set)
        for goid_child, goid_parent in self.edges:
            c2ps[goid_child].add(goid_parent)
        return c2ps

    def _getobjs_higher(self, goobj):
        """Get all parents/relationships on this GOTerm."""
        goobjs_higher = set(goobj.parents)
        for reltyp, relgoobjs in goobj.relationship.items():
            if reltyp in self.relationships:
                goobjs_higher.update(relgoobjs)
        return goobjs_higher



# -- Initialization by considering all child and/or parent relatives -----------
class EdgesRelatives(EdgesBase):
    """Inits GO-to-GO edges using all relatives above and/or below source GOs."""

    # pylint: disable=too-many-arguments
    # def __init__(self, go2obj, relationships, go_sources, traverse_parent, traverse_child):
    def __init__(self, gosubdag, traverse_parent, traverse_child):
        super(EdgesRelatives, self).__init__(gosubdag)
        # go2obj contain GO IDs in subset
        _gos = set(gosubdag.go2obj)
        assert traverse_child or traverse_parent, "NO EDGES IN GRAPH"
        # GO IDs for child->parents
        p2cs = self._init_p2cs(_gos, traverse_parent)
        # GO IDs for parent->children
        c2ps = self._init_c2ps(gosubdag.go_sources, traverse_child)
        # GO IDs for GO->relationship
        rel2src2dsts = self._init_rel2src2dsts(_gos, traverse_parent)
        rel2dst2srcs = self._init_rel2dst2srcs(_gos, traverse_child)
        # Set by derived edge class
        # self.edges = self._init_edges(_gos, p2cs, c2ps)
        self.edges = self._init_edges(p2cs, c2ps)
        self.edges_rel = self._init_edges_relationships(rel2src2dsts, rel2dst2srcs)
        assert _gos == set(self.go2obj)
        # self.chk_edges()

    @staticmethod
    # Too slow to check goids_present as we go. Only minor init modes need checking.
    # def _init_edges(goids_present, p2cs, c2ps):
    def _init_edges(p2cs, c2ps):
        """Get the directed edges from GO term to GO term."""
        edge_from_to = []
        for parent, children in p2cs.items():
            for child in children:
                # if child in goids_present and parent in goids_present:
                edge_from_to.append((child, parent))
        for parent, children in c2ps.items():
            for child in children:
                # if child in goids_present and parent in goids_present:
                edge_from_to.append((child, parent))
        return edge_from_to

    @staticmethod
    def _init_edges_relationships(rel2src2dsts, rel2dst2srcs):
        """Get the directed edges from GO term to GO term using relationships."""
        edge_rel2fromto = {}
        relationships = set(rel2src2dsts).union(rel2dst2srcs)
        for reltype in relationships:
            edge_from_to = []
            if reltype in rel2src2dsts:
                for parent, children in rel2src2dsts[reltype].items():
                    for child in children:
                        edge_from_to.append((child, parent))
            if reltype in rel2dst2srcs:
                for parent, children in rel2dst2srcs[reltype].items():
                    for child in children:
                        edge_from_to.append((child, parent))
            edge_rel2fromto[reltype] = edge_from_to
        return edge_rel2fromto

    # -------------------------------------------------------------------
    def _init_rel2src2dsts(self, go_sources, traverse_parent):
        """Traverse up parents."""
        if not traverse_parent or not self.relationships:
            return {}
        rel2src2dsts = {r:defaultdict(set) for r in self.relationships}
        goids_seen = set()
        go2obj = self.go2obj
        for goid_src in go_sources:
            goobj_src = go2obj[goid_src]
            if goobj_src.relationship and goid_src not in goids_seen:
                self._traverse_relationship_objs(rel2src2dsts, goobj_src, goids_seen)
        return rel2src2dsts

    def _traverse_relationship_objs(self, rel2src2dsts, goobj_child, goids_seen):
        """Traverse from source GO up relationships."""
        child_id = goobj_child.id
        goids_seen.add(child_id)
        ##A self.go2obj[child_id] = goobj_child
        # Update goids_seen and go2obj with child alt_ids
        for goid_altid in goobj_child.alt_ids:
            goids_seen.add(goid_altid)
            ##A self.go2obj[goid_altid] = goobj_child
        # Loop through relationships of child object
        for reltype, recs in goobj_child.relationship.items():
            if reltype in self.relationships:
                for relationship_obj in recs:
                    relationship_id = relationship_obj.id
                    rel2src2dsts[reltype][relationship_id].add(child_id)
                    # If relationship has not been seen, traverse
                    if relationship_id not in goids_seen:
                        self._traverse_relationship_objs(rel2src2dsts, relationship_obj, goids_seen)

    # -------------------------------------------------------------------
    def _init_rel2dst2srcs(self, go_sources, traverse_child):
        """Traverse through reverse relationships."""
        if not traverse_child or not self.relationships:
            return {}
        rel2dst2srcs = {r:defaultdict(set) for r in self.relationships}
        goids_seen = set()
        go2obj = self.go2obj
        for goid_src in go_sources:
            goobj_src = go2obj[goid_src]
            if goid_src not in goids_seen:
                self._traverse_relationship_rev_objs(rel2dst2srcs, goobj_src, goids_seen)
        return rel2dst2srcs

    def _traverse_relationship_rev_objs(self, rel2dst2srcs, goobj_parent, goids_seen):
        """Traverse from source GO down children."""
        parent_id = goobj_parent.id
        goids_seen.add(parent_id)
        ##A self.go2obj[parent_id] = goobj_parent
        # Update goids_seen and go2obj with parent alt_ids
        for goid_altid in goobj_parent.alt_ids:
            goids_seen.add(goid_altid)
            ##A self.go2obj[goid_altid] = goobj_parent
        # Loop through children
        for reltype, recs in goobj_parent.relationship.items():
            if reltype in self.relationships:
                for relrev_obj in recs:
                    relrev_id = relrev_obj.id
                    rel2dst2srcs[relrev_id].add(parent_id)
                    # If child has not been seen, traverse
                    if relrev_id not in goids_seen:
                        ##F self._traverse_relrev_objs(rel2dst2srcs, relrev_obj, go2obj, goids_seen)
                        self._traverse_relationship_rev_objs(rel2dst2srcs, relrev_obj, goids_seen)

    # -------------------------------------------------------------------
    def _init_p2cs(self, go_sources, traverse_parent):
        """Traverse up parents."""
        if not traverse_parent:
            return {}
        p2cs = defaultdict(set)
        goids_seen = set()
        go2obj = self.go2obj
        for goid_src in go_sources:
            goobj_src = go2obj[goid_src]
            if goid_src not in goids_seen:
                ##F self._traverse_parent_objs(p2cs, goobj_src, go2obj, goids_seen)
                self._traverse_parent_objs(p2cs, goobj_src, goids_seen)
        return p2cs

    ##F def _traverse_parent_objs(self, p2cs, goobj_child, go2obj, goids_seen):
    def _traverse_parent_objs(self, p2cs, goobj_child, goids_seen):
        """Traverse from source GO up parents."""
        # Update public(go2obj p2cs), private(goids_seen)
        child_id = goobj_child.id
        # mark child as seen
        goids_seen.add(child_id)
        ##A self.go2obj[child_id] = goobj_child
        # Update goids_seen and go2obj with child alt_ids
        for goid_altid in goobj_child.alt_ids:
            goids_seen.add(goid_altid)
            ##A self.go2obj[goid_altid] = goobj_child
        # Loop through parents of child object
        for parent_obj in goobj_child.parents:
            parent_id = parent_obj.id
            p2cs[parent_id].add(child_id)
            # If parent has not been seen, traverse
            if parent_id not in goids_seen:
                ##F self._traverse_parent_objs(p2cs, parent_obj, go2obj, goids_seen)
                self._traverse_parent_objs(p2cs, parent_obj, goids_seen)

    # -------------------------------------------------------------------
    def _init_c2ps(self, go_sources, traverse_child):
        """Traverse up children."""
        if not traverse_child:
            return {}
        c2ps = defaultdict(set)
        goids_seen = set()
        go2obj = self.go2obj
        for goid_src in go_sources:
            goobj_src = go2obj[goid_src]
            if goid_src not in goids_seen:
                ##F self._traverse_child_objs(c2ps, goobj_src, go2obj, goids_seen)
                self._traverse_child_objs(c2ps, goobj_src, goids_seen)
        return c2ps

    ##F def _traverse_child_objs(self, c2ps, goobj_parent, go2obj, goids_seen):
    def _traverse_child_objs(self, c2ps, goobj_parent, goids_seen):
        """Traverse from source GO down children."""
        # Update public(godag.go2obj godag.c2ps), private(_seen_pids)
        parent_id = goobj_parent.id
        # mark parent as seen
        goids_seen.add(parent_id)
        ##A self.go2obj[parent_id] = goobj_parent
        # Update goids_seen and go2obj with parent alt_ids
        for goid_altid in goobj_parent.alt_ids:
            goids_seen.add(goid_altid)
            ##A self.go2obj[goid_altid] = goobj_parent
        # Loop through children
        for child_obj in goobj_parent.children:
            child_id = child_obj.id
            c2ps[child_id].add(parent_id)
            # If child has not been seen, traverse
            if child_id not in goids_seen:
                ##F self._traverse_child_objs(c2ps, child_obj, go2obj, goids_seen)
                self._traverse_child_objs(c2ps, child_obj, goids_seen)


# -- Initialization with realtives on specific src-dst paths -------------------
class EdgesPath(EdgesBase):
    """Inits GO-to-GO edges using a list of (parent destination, child sources)"""

    def __init__(self, gosubdag, dst_srcs_list):
        super(EdgesPath, self).__init__(gosubdag)
        self.edges = None
        self.goid_all = None
        self._init_edges(dst_srcs_list)
        # GO IDs for child->parents
        # self.p2cs = self._init_p2cs(go_sources, traverse_parent)
        # GO IDs for parent->children
        # self.c2ps = self._init_c2ps(go_sources, traverse_child)

    def get_edges(self):
        """Get the directed edges from GO term to GO term."""
        return self.edges

    def _init_edges(self, dst_srcs_list):
        """Create all GO edges given a list of (dst, srcs)."""
        from goatools.gosubdag.go_paths import get_paths_goobjs, paths2edges
        edges_all = set()
        goid_all = set()
        go2obj = self.go2obj
        for dst, srcs in dst_srcs_list:
            go2obj_srcs = {}
            for goid in srcs:
                go2obj_srcs[goid] = go2obj[goid]
            go_paths, go_all = get_paths_goobjs(go2obj_srcs.values(), go_top=dst, go2obj=go2obj)
            edges_all |= paths2edges(go_paths)
            goid_all |= go_all
        self.edges = [(a.id, b.id) for a, b in edges_all]
        self.goid_all = goid_all

# Copyright (C) 2016-2018, DV Klopfenstein, H Tang, All rights reserved.
