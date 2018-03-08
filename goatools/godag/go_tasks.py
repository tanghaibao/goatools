"""GO-DAG tasks."""

__copyright__ = "Copyright (C) 2010-2018, DV Klopfenstein, H Tang, All rights reserved."
__author__ = "DV Klopfenstein"

# ------------------------------------------------------------------------------------
def get_go2parents(goterms):
    """Get all parent GO IDs for each GO in dict keys."""
    go2parents = {}
    for goterm in goterms:
        _get_go2parents(go2parents, goterm.id, goterm)
    return go2parents

def get_go2children(goterms):
    """Get all parent GO IDs for each GO in dict keys."""
    go2children = {}
    for goterm in goterms:
        _get_go2children(go2children, goterm.id, goterm)
    return go2children

def get_go2upper(goterms):
    """Get all parent GO IDs for each GO in dict keys."""
    go2upper = {}
    for goterm in goterms:
        _get_go2upper(go2upper, goterm.id, goterm)
    return go2upper

def get_go2lower(goterms):
    """Get all parent GO IDs for each GO in dict keys."""
    go2lower = {}
    for goterm in goterms:
        _get_go2lower(go2lower, goterm.id, goterm)
    return go2lower

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
def _get_go2parents(go2parents, goid, goterm):
    """Add the parent GO IDs for one GO term and their parents."""
    if goid in go2parents:
        return go2parents[goid]
    parent_goids = set()
    for parent_goterm in goterm.parents:
        parent_goid = parent_goterm.id
        parent_goids.add(parent_goid)
        parent_goids |= _get_go2parents(go2parents, parent_goid, parent_goterm)
    go2parents[goid] = parent_goids
    return parent_goids

def _get_go2children(go2children, goid, goterm):
    """Add the child GO IDs for one GO term and their children."""
    if goid in go2children:
        return go2children[goid]
    child_goids = set()
    for child_goterm in goterm.children:
        child_goid = child_goterm.id
        child_goids.add(child_goid)
        child_goids |= _get_go2children(go2children, child_goid, child_goterm)
    go2children[goid] = child_goids
    return child_goids

def _get_go2upper(go2upper, goid, goterm):
    """Add the parent GO IDs for one GO term and their upper."""
    if goid in go2upper:
        return go2upper[goid]
    upper_goids = set()
    for upper_goterm in goterm.get_goterms_upper():
        upper_goid = upper_goterm.id
        upper_goids.add(upper_goid)
        upper_goids |= _get_go2upper(go2upper, upper_goid, upper_goterm)
    go2upper[goid] = upper_goids
    return upper_goids

def _get_go2lower(go2lower, goid, goterm):
    """Add the lower GO IDs for one GO term and their lowerren."""
    if goid in go2lower:
        return go2lower[goid]
    lower_goids = set()
    for lower_goterm in goterm.get_goterms_lower():
        lower_goid = lower_goterm.id
        lower_goids.add(lower_goid)
        lower_goids |= _get_go2lower(go2lower, lower_goid, lower_goterm)
    go2lower[goid] = lower_goids
    return lower_goids

# ------------------------------------------------------------------------------------
class CurNHigher(object):
    """Fill go2obj with GO IDs in relationships."""

    def __init__(self, relationships, go2obj_all):
        # True or A set of relationships we would like to keep
        self.relationships = relationships
        self.go2obj_all = go2obj_all

    def get_go2obj_cur_n_high(self, go2obj_user, go_sources):
        """Get go2obj containing: go_srcs and parents."""
        if not self.relationships:
            self._get_go2obj_high(go2obj_user, go_sources, self.fill_parentgoid2obj_r0)
        else:
            self._get_go2obj_high(go2obj_user, go_sources, self.fill_parentgoid2obj_r1)

    def _get_go2obj_high(self, go2obj_user, go_sources, fnc_fill):
        """Get go2obj containing: go_srcs and parents."""
        for goid_user in go_sources:
            goobj_user = self.go2obj_all[goid_user]
            fnc_fill(go2obj_user, goobj_user)
            go2obj_user[goobj_user.id] = goobj_user
            if goid_user != goobj_user.id:
                go2obj_user[goid_user] = goobj_user

    def fill_parentgoid2obj_r0(self, go2obj, child_obj):
        """Fill go2obj with all parent key GO IDs and their objects."""
        for parent_obj in child_obj.parents:
            if parent_obj.id not in go2obj:
                go2obj[parent_obj.id] = parent_obj
                self.fill_parentgoid2obj_r0(go2obj, parent_obj)

    def fill_parentgoid2obj_r1(self, go2obj_user, child_obj):
        """Fill go2obj_user with all parent/relationship key GO IDs and their objects."""
        for higher_obj in self._getobjs_higher(child_obj):
            if higher_obj.id not in go2obj_user:
                go2obj_user[higher_obj.id] = higher_obj
                self.fill_parentgoid2obj_r1(go2obj_user, higher_obj)

    def _getobjs_higher(self, goobj):
        """Get all parents/relationships on this GOTerm."""
        goobjs_higher = set(goobj.parents)
        for reltyp, relgoobjs in goobj.relationship.items():
            if self.relationships is True or reltyp in self.relationships:
                goobjs_higher.update(relgoobjs)
        return goobjs_higher


# Copyright (C) 2010-2018, DV Klopfenstein, H Tang, All rights reserved.
