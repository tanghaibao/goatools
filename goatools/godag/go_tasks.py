"""GO-DAG tasks."""

__copyright__ = "Copyright (C) 2010-2018, DV Klopfenstein, H Tang, All rights reserved."
__author__ = "DV Klopfenstein"

# ------------------------------------------------------------------------------------
def get_go2parents(goterms):
    """Get all parent GO IDs for each GO in dict keys."""
    go2parents = {}
    for goterm in goterms:
        _get_go2parents(goterm.id, goterm, go2parents)
    return go2parents

def get_go2children(goterms):
    """Get all parent GO IDs for each GO in dict keys."""
    go2children = {}
    for goterm in goterms:
        _get_go2children(goterm.id, goterm, go2children)
    return go2children

def get_relationship_targets(goids, relationships, go2rec):
    """Get GO ID set of GO IDs in a relationship target set."""
    # Requirements to use this function:
    #     1) GO Terms must have been loaded with 'relationships'
    #     2) GO IDs in 'goids' arguement must be present in go2rec
    #     3) Arg, 'relationships' must be True or an iterable
    reltgt_goterms_all = set()
    for goid in goids:
        goterm = go2rec[goid]
        for reltype, reltgt_goterms_cur in goterm.relationship.items():
            if relationships is True or reltype in relationships:
                reltgt_goterms_all.update(reltgt_goterms_cur)
    return reltgt_goterms_all

# ------------------------------------------------------------------------------------
def _get_go2parents(goid, goterm, go2parents):
    """Add the parent GO IDs for one GO term and their parents."""
    if goid in go2parents:
        return go2parents[goid]
    parent_goids = set()
    for parent_goterm in goterm.parents:
        parent_goid = parent_goterm.id
        parent_goids.add(parent_goid)
        parent_goids |= _get_go2parents(parent_goid, parent_goterm, go2parents)
    go2parents[goid] = parent_goids
    return parent_goids

def _get_go2children(goid, goterm, go2children):
    """Add the child GO IDs for one GO term and their children."""
    if goid in go2children:
        return go2children[goid]
    child_goids = set()
    for child_goterm in goterm.children:
        child_goid = child_goterm.id
        child_goids.add(child_goid)
        child_goids |= _get_go2children(child_goid, child_goterm, go2children)
    go2children[goid] = child_goids
    return child_goids


# Copyright (C) 2010-2018, DV Klopfenstein, H Tang, All rights reserved.
