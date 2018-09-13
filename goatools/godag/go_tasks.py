"""item-DAG tasks."""

__copyright__ = "Copyright (C) 2010-2018, DV Klopfenstein, H Tang, All rights reserved."
__author__ = "DV Klopfenstein"

# ------------------------------------------------------------------------------------
def get_id2parents(objs):
    """Get all parent item IDs for each item in dict keys."""
    id2parents = {}
    for obj in objs:
        _get_id2parents(id2parents, obj.item_id, obj)
    return id2parents

def get_id2children(objs):
    """Get all parent item IDs for each item in dict keys."""
    id2children = {}
    for obj in objs:
        _get_id2children(id2children, obj.item_id, obj)
    return id2children

def get_id2upper(objs):
    """Get all parent item IDs for each item in dict keys."""
    id2upper = {}
    for obj in objs:
        _get_id2upper(id2upper, obj.item_id, obj)
    return id2upper

def get_id2lower(objs):
    """Get all parent item IDs for each item in dict keys."""
    id2lower = {}
    for obj in objs:
        _get_id2lower(id2lower, obj.item_id, obj)
    return id2lower

def get_relationship_targets(item_ids, relationships, id2rec):
    """Get item ID set of item IDs in a relationship target set."""
    # Requirements to use this function:
    #     1) item Terms must have been loaded with 'relationships'
    #     2) item IDs in 'item_ids' arguement must be present in id2rec
    #     3) Arg, 'relationships' must be True or an iterable
    reltgt_objs_all = set()
    for goid in item_ids:
        obj = id2rec[goid]
        for reltype, reltgt_objs_cur in obj.relationship.items():
            if relationships is True or reltype in relationships:
                reltgt_objs_all.update(reltgt_objs_cur)
    return reltgt_objs_all

# ------------------------------------------------------------------------------------
def _get_id2parents(id2parents, item_id, item_obj):
    """Add the parent item IDs for one item object and their parents."""
    if item_id in id2parents:
        return id2parents[item_id]
    parent_ids = set()
    for parent_obj in item_obj.parents:
        parent_id = parent_obj.item_id
        parent_ids.add(parent_id)
        parent_ids |= _get_id2parents(id2parents, parent_id, parent_obj)
    id2parents[item_id] = parent_ids
    return parent_ids

def _get_id2children(id2children, item_id, item_obj):
    """Add the child item IDs for one item object and their children."""
    if item_id in id2children:
        return id2children[item_id]
    child_ids = set()
    for child_obj in item_obj.children:
        child_id = child_obj.item_id
        child_ids.add(child_id)
        child_ids |= _get_id2children(id2children, child_id, child_obj)
    id2children[item_id] = child_ids
    return child_ids

def _get_id2upper(id2upper, item_id, item_obj):
    """Add the parent item IDs for one item object and their upper."""
    if item_id in id2upper:
        return id2upper[item_id]
    upper_ids = set()
    for upper_obj in item_obj.get_goterms_upper():
        upper_id = upper_obj.item_id
        upper_ids.add(upper_id)
        upper_ids |= _get_id2upper(id2upper, upper_id, upper_obj)
    id2upper[item_id] = upper_ids
    return upper_ids

def _get_id2lower(id2lower, item_id, item_obj):
    """Add the lower item IDs for one item object and the objects below them."""
    if item_id in id2lower:
        return id2lower[item_id]
    lower_ids = set()
    for lower_obj in item_obj.get_goterms_lower():
        lower_id = lower_obj.item_id
        lower_ids.add(lower_id)
        lower_ids |= _get_id2lower(id2lower, lower_id, lower_obj)
    id2lower[item_id] = lower_ids
    return lower_ids

# ------------------------------------------------------------------------------------
class CurNHigher(object):
    """Fill id2obj with item IDs in relationships."""

    def __init__(self, relationships, id2obj_all):
        # True or A set of relationships we would like to keep
        self.relationships = relationships
        self.id2obj_all = id2obj_all

    def get_id2obj_cur_n_high(self, id2obj_user, id_sources):
        """Get id2obj containing: id_srcs and parents."""
        if not self.relationships:
            self._get_id2obj_high(id2obj_user, id_sources, self.fill_parentidid2obj_r0)
        else:
            self._get_id2obj_high(id2obj_user, id_sources, self.fill_parentidid2obj_r1)

    def _get_id2obj_high(self, id2obj_user, id_sources, fnc_fill):
        """Get id2obj containing: id_srcs and parents."""
        for idid_user in id_sources:
            idobj_user = self.id2obj_all[idid_user]
            fnc_fill(id2obj_user, idobj_user)
            id2obj_user[idobj_user.item_id] = idobj_user
            if idid_user != idobj_user.item_id:
                id2obj_user[idid_user] = idobj_user

    def fill_parentidid2obj_r0(self, id2obj, child_obj):
        """Fill id2obj with all parent key item IDs and their objects."""
        for parent_obj in child_obj.parents:
            if parent_obj.item_id not in id2obj:
                id2obj[parent_obj.item_id] = parent_obj
                self.fill_parentidid2obj_r0(id2obj, parent_obj)

    def fill_parentidid2obj_r1(self, id2obj_user, child_obj):
        """Fill id2obj_user with all parent/relationship key item IDs and their objects."""
        for higher_obj in self._getobjs_higher(child_obj):
            if higher_obj.item_id not in id2obj_user:
                id2obj_user[higher_obj.item_id] = higher_obj
                self.fill_parentidid2obj_r1(id2obj_user, higher_obj)

    def _getobjs_higher(self, idobj):
        """Get all parents/relationships on this GOTerm."""
        idobjs_higher = set(idobj.parents)
        for reltyp, relidobjs in idobj.relationship.items():
            if self.relationships is True or reltyp in self.relationships:
                idobjs_higher.update(relidobjs)
        return idobjs_higher


# Copyright (C) 2010-2018, DV Klopfenstein, H Tang, All rights reserved.
