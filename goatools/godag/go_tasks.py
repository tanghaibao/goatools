"""GO-DAG tasks."""

__copyright__ = "Copyright (C) 2010-2018, DV Klopfenstein, H Tang, All rights reserved."
__author__ = "DV Klopfenstein"

# ------------------------------------------------------------------------------------
def get_go2parents(goobjs):
    """Given go2obj, get all parent GO IDs for each GO in dict keys."""
    go2parents = {}
    for goobj in goobjs:
        _get_go2parents(goobj.id, goobj, go2parents)
    return go2parents

def get_go2children(goobjs):
    """Given go2obj, get all parent GO IDs for each GO in dict keys."""
    go2children = {}
    for goobj in goobjs:
        _get_go2children(goobj.id, goobj, go2children)
    return go2children


# ------------------------------------------------------------------------------------
def _get_go2parents(goid, goobj, go2parents):
    """Add the parent GO IDs for one GO term and their parents."""
    if goid in go2parents:
        return go2parents[goid]
    parent_goids = set()
    for parent_goobj in goobj.parents:
        parent_goid = parent_goobj.id
        parent_goids.add(parent_goid)
        parent_goids |= _get_go2parents(parent_goid, parent_goobj, go2parents)
    go2parents[goid] = parent_goids
    return parent_goids

def _get_go2children(goid, goobj, go2children):
    """Add the child GO IDs for one GO term and their children."""
    if goid in go2children:
        return go2children[goid]
    child_goids = set()
    for child_goobj in goobj.children:
        child_goid = child_goobj.id
        child_goids.add(child_goid)
        child_goids |= _get_go2children(child_goid, child_goobj, go2children)
    go2children[goid] = child_goids
    return child_goids


# Copyright (C) 2010-2018, DV Klopfenstein, H Tang, All rights reserved.
