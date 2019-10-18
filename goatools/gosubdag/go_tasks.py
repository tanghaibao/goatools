"""Tasks for go2obj dicts."""

__copyright__ = "Copyright (C) 2016-present, DV Klopfenstein, H Tang, All rights reserved."
__author__ = "DV Klopfenstein"

import collections as cx
from goatools.godag.go_tasks import get_go2ancestors
from goatools.godag.go_tasks import get_go2descendants


# ------------------------------------------------------------------------------------
def get_sorted_relationship(goterms):
    """Topological sort of GO Terms w/'relationship's loaded."""
    return TopologicalSortRelationships(goterms).goterms_sorted

class TopologicalSortRelationships:
    """Topological sort of GO Terms w/'relationship's loaded."""

    # pylint: disable=too-few-public-methods
    def __init__(self, goterms):
        self.goterms_sorted = []
        self.goids_seen = set()
        self._init_sorted_relationship(goterms)

    def _init_sorted_relationship(self, goterms):
        """Topologically sort GO Terms using 'is_a' parents and 'relationship' GO IDs."""
        # NOTE: GODag must be loaded with 'relationship' to use this function
        for goterm in goterms:
            self._get_sorted_relationships(goterm)

    def _get_sorted_relationships(self, goterm):
        """Traverse GO Terms above the current GO Term. Then add current GO Term to sorted."""
        if goterm.id in self.goids_seen:
            return
        self.goids_seen.add(goterm.id)
        for goterm_upper in goterm.get_goterms_upper():
            self._get_sorted_relationships(goterm_upper)
        self.goterms_sorted.append(goterm)

# ------------------------------------------------------------------------------------
def get_go2obj_unique(go2obj):
    """If GO keys point to the same GOTerm, return new go2obj w/no duplicates."""
    # Find the unique GO Terms that are represented for each GO in go2obj
    goid2gokeys = cx.defaultdict(set)
    for goid, goobj in go2obj.items():
        goid2gokeys[goobj.id].add(goid)
    go_unique = set()
    for goid, gos_seen in goid2gokeys.items():
        # Return main GO ID, if it is present in the go2obj keys
        if goid in gos_seen:
            go_unique.add(goid)
        # Otherwise return an alternate GO ID
        else:
            go_unique.add(next(iter(gos_seen)))
    return go_unique

# ------------------------------------------------------------------------------------
def get_go2parents_go2obj(go2obj, relationships=None, prt=None):
    """Return go2ancestors (set of parent GO IDs) for all GO ID keys in go2obj."""
    goobjs, altgo2goobj = get_goobjs_altgo2goobj(go2obj)
    go2ancestors = get_go2ancestors(goobjs, relationships, prt)
    add_alt_goids(go2ancestors, altgo2goobj)
    return go2ancestors

# ------------------------------------------------------------------------------------
def get_go2children_go2obj(go2obj, relationships=None, prt=None):
    """Return go2children (set of child GO IDs) for all GO ID keys in go2obj."""
    goobjs, altgo2goobj = get_goobjs_altgo2goobj(go2obj)
    go2children = get_go2descendants(goobjs, relationships, prt)
    add_alt_goids(go2children, altgo2goobj)
    return go2children

# ------------------------------------------------------------------------------------
def get_goobjs_altgo2goobj(go2obj):
    """Separate alt GO IDs and key GO IDs."""
    goobjs = set()
    altgo2goobj = {}
    for goid, goobj in go2obj.items():
        goobjs.add(goobj)
        if goid != goobj.id:
            altgo2goobj[goid] = goobj
    return goobjs, altgo2goobj

def add_alt_goids(go2values, altgo2goobj):
    """Add alternate source GO IDs."""
    for goobj_key in altgo2goobj.values():
        if goobj_key.id in go2values:
            values_curr = go2values[goobj_key.id]
            for goid_alt in goobj_key.alt_ids:
                go2values[goid_alt] = values_curr
    return go2values

# ------------------------------------------------------------------------------------
def fill_main_goids(go2obj, goids):
    """Ensure main GO IDs are included in go2obj."""
    # User GO IDs (goids) may be either main GO IDs or alternate GO IDs.
    for goid in goids:
        goobj = go2obj[goid]
        # If a user specified an ALT GO ID and main GO ID not in go2obj:
        if goid != goobj.id and goobj.id not in go2obj:
            # Add main GO ID to go2obj
            go2obj[goobj.id] = goobj

def fill_altgoids(go2obj):
    """Given a go2obj containing key GO IDs, fill with all alternate GO IDs."""
    alt2obj = {altgo:goobj for goobj in go2obj.values() for altgo in goobj.alt_ids}
    for goid, goobj in alt2obj.items():
        go2obj[goid] = goobj

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
def fill_relationshipobjs(go2obj, relationships):
    """Add GO IDs to go2obj that are involved in relationships."""
    # Get all GO Term record objects that have relationships
    obj = RelationshipFill(go2obj, relationships)
    for goobj in go2obj.values():
        if goobj.relationship:
            obj.fill_relationshipgo2obj(goobj)
        if goobj.relationship_rev:
            obj.fill_relationshiprevgo2obj(goobj)

class RelationshipFill:
    """Fill go2obj with GO IDs in relatinships."""

    def __init__(self, go2obj, relationships):
        # This dict shall be augmented with higher parent/relationship GO IDs
        self.go2obj = go2obj
        # A set of relationships we would like to keep
        self.relationships = relationships

    def fill_relationshipgo2obj(self, goobj):
        """Fill go2obj with all relationship key GO IDs and their objects."""
        for reltyp, relgoobjs in goobj.relationship.items():
            if reltyp in self.relationships:
                for relgoobj in relgoobjs:
                    if relgoobj.id not in self.go2obj:
                        self.go2obj[relgoobj.id] = relgoobj
                        self.fill_relationshipgo2obj(relgoobj)

    def fill_relationshiprevgo2obj(self, goobj):
        """Fill go2obj with all relationship key GO IDs and their objects."""
        for reltyp, relgoobjs in goobj.relationship_rev.items():
            if reltyp in self.relationships:
                for relgoobj in relgoobjs:
                    if relgoobj.id not in self.go2obj:
                        self.go2obj[relgoobj.id] = relgoobj
                        self.fill_relationshiprevgo2obj(relgoobj)

# ------------------------------------------------------------------------------------
def get_child_objs(parent_obj):
    """Fill child2obj with all child key and alt GO IDs and their objects."""
    child2obj = {}
    fill_childgoid2obj(child2obj, parent_obj)
    fill_altgoids(child2obj)
    return child2obj

def fill_childgoid2obj(childgoid2obj, parent_obj):
    """Fill childgoid2obj with all child key GO IDs and their objects."""
    for child_obj in parent_obj.children:
        if child_obj.id not in childgoid2obj:
            childgoid2obj[child_obj.id] = child_obj
            fill_childgoid2obj(childgoid2obj, child_obj)

# ------------------------------------------------------------------------------------
def get_leaf_children(gos_user, go2obj_arg):
    """Find all the GO descendants under all user GO IDs. Return leaf-level GO IDs."""
    childgoid2obj = {}
    for goid_usr in gos_user:
        goobj_usr = go2obj_arg[goid_usr]
        fill_childgoid2obj(childgoid2obj, goobj_usr)
    return set(go for go, o in childgoid2obj.items() if not o.children)

# ------------------------------------------------------------------------------------
def goid_is_valid(goid):
    """Check format of user-provided GO IDs"""
    return goid[:3] == "GO:" and len(goid) == 10 and goid[3:].isdigit()

def goids_valid(goids):
    """Check format of user-provided GO IDs"""
    for goid in goids:
        if not goid_is_valid(goid):
            return False
    return True

def chk_goids(goids, msg=None, raise_except=True):
    """check that all GO IDs have the proper format."""
    for goid in goids:
        if not goid_is_valid(goid):
            if raise_except:
                raise RuntimeError("BAD GO({GO}): {MSG}".format(GO=goid, MSG=msg))
            return goid
    return None


# Copyright (C) 2016-present, DV Klopfenstein, H Tang, All rights reserved.
