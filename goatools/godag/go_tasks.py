"""item-DAG tasks."""

__copyright__ = "Copyright (C) 2010-present, DV Klopfenstein, H Tang, All rights reserved."
__author__ = "DV Klopfenstein"

from goatools.godag.consts import RELATIONSHIP_SET


# ------------------------------------------------------------------------------------
def get_go2parents(go2obj, relationships):
    """Get set of parents GO IDs, including parents through user-specfied relationships"""
    if go2obj and not hasattr(next(iter(go2obj.values())), 'relationship') or not relationships:
        return get_go2parents_isa(go2obj)
    go2parents = {}
    for goid_main, goterm in go2obj.items():
        parents_goids = set(o.id for o in goterm.parents)
        for rel in set(goterm.relationship.keys()).intersection(relationships):
            for parent_term in goterm.relationship[rel]:
                parents_goids.add(parent_term.id)
        if parents_goids:
            go2parents[goid_main] = parents_goids
    return go2parents

# ------------------------------------------------------------------------------------
def get_go2children(go2obj, relationships):
    """Get set of children GO IDs, including children through user-specfied relationships"""
    if go2obj and not hasattr(next(iter(go2obj.values())), 'relationship') or not relationships:
        return get_go2children_isa(go2obj)
    go2children = {}
    for goid_main, goterm in go2obj.items():
        children_goids = set(o.id for o in goterm.children)
        for rel in set(goterm.relationship_rev.keys()).intersection(relationships):
            for child_term in goterm.relationship_rev[rel]:
                children_goids.add(child_term.id)
        if children_goids:
            go2children[goid_main] = children_goids
    return go2children

# ------------------------------------------------------------------------------------
def get_go2parents_isa(go2obj):
    """Get set of immediate parents GO IDs"""
    go2parents = {}
    for goid_main, goterm in go2obj.items():
        parents_goids = set(o.id for o in goterm.parents)
        if parents_goids:
            go2parents[goid_main] = parents_goids
    return go2parents

# ------------------------------------------------------------------------------------
def get_go2children_isa(go2obj):
    """Get set of immediate children GO IDs"""
    go2children = {}
    for goid_main, goterm in go2obj.items():
        children_goids = set(o.id for o in goterm.children)
        if children_goids:
            go2children[goid_main] = children_goids
    return go2children

# ------------------------------------------------------------------------------------
def get_go2ancestors(terms, relationships, prt=None):
    """Get GO-to- ancestors (all parents)"""
    if not relationships:
        if prt is not None:
            prt.write('up: is_a\n')
        return get_id2parents(terms)
    if relationships == RELATIONSHIP_SET or relationships is True:
        if prt is not None:
            prt.write('up: is_a and {Rs}\n'.format(
                Rs=' '.join(sorted(RELATIONSHIP_SET))))
        return get_id2upper(terms)
    if prt is not None:
        prt.write('up: is_a and {Rs}\n'.format(
            Rs=' '.join(sorted(relationships))))
    return get_id2upperselect(terms, relationships)

def get_go2descendants(terms, relationships, prt=None):
    """Get GO-to- descendants"""
    if not relationships:
        if prt is not None:
            prt.write('down: is_a\n')
        return get_id2children(terms)
    if relationships == RELATIONSHIP_SET or relationships is True:
        if prt is not None:
            prt.write('down: is_a and {Rs}\n'.format(
                Rs=' '.join(sorted(RELATIONSHIP_SET))))
        return get_id2lower(terms)
    if prt is not None:
        prt.write('down: is_a and {Rs}\n'.format(
            Rs=' '.join(sorted(relationships))))
    return get_id2lowerselect(terms, relationships)

# ------------------------------------------------------------------------------------
def get_id2parents(objs):
    """Get all parent IDs up the hierarchy"""
    id2parents = {}
    for obj in objs:
        _get_id2parents(id2parents, obj.item_id, obj)
    return {e:es for e, es in id2parents.items() if es}

def get_id2children(objs):
    """Get all child IDs down the hierarchy"""
    id2children = {}
    for obj in objs:
        _get_id2children(id2children, obj.item_id, obj)
    return {e:es for e, es in id2children.items() if es}

def get_id2upper(objs):
    """Get all ancestor IDs, including all parents and IDs up all relationships"""
    id2upper = {}
    for obj in objs:
        _get_id2upper(id2upper, obj.item_id, obj)
    return {e:es for e, es in id2upper.items() if es}

def get_id2lower(objs):
    """Get all descendant IDs, including all children and IDs down all relationships"""
    id2lower = {}
    for obj in objs:
        _get_id2lower(id2lower, obj.item_id, obj)
    return {e:es for e, es in id2lower.items() if es}

def get_id2upperselect(objs, relationship_set):
    """Get all ancestor IDs, including all parents and IDs up selected relationships"""
    return IdToUpperSelect(objs, relationship_set).id2upperselect

def get_id2lowerselect(objs, relationship_set):
    """Get all descendant IDs, including all children and IDs down selected relationships"""
    return IdToLowerSelect(objs, relationship_set).id2lowerselect

def get_relationship_targets(item_ids, relationships, id2rec):
    """Get item ID set of item IDs in a relationship target set"""
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
# pylint: disable=too-few-public-methods
class IdToUpperSelect:
    """Get all ancestor IDs, including all parents and IDs up selected relationships"""

    def __init__(self, objs, relationship_set):
        self.rset = relationship_set
        self.id2upperselect = {}
        self._init_id2upperselect(objs)

    def _init_id2upperselect(self, objs):
        """Get all parent item IDs for each item in dict keys."""
        for obj in objs:
            self._get_id2upperselect(obj.item_id, obj)

    def _get_id2upperselect(self, item_id, item_obj):
        """Add the parent item IDs for one item object and their parents."""
        id2upperselect = self.id2upperselect
        if item_id in id2upperselect:
            return id2upperselect[item_id]
        parent_ids = set()
        r2os = item_obj.relationship
        sets = [r2os[r] for r in self.rset.intersection(r2os)]
        for parent_obj in item_obj.parents.union(*sets):
            parent_id = parent_obj.item_id
            parent_ids.add(parent_id)
            parent_ids |= self._get_id2upperselect(parent_id, parent_obj)
        id2upperselect[item_id] = parent_ids
        return parent_ids

class IdToLowerSelect:
    """Get all descendant IDs, including all children and IDs down selected relationships"""

    def __init__(self, objs, relationship_set):
        self.rset = relationship_set
        self.id2lowerselect = {}
        self._init_id2lowerselect(objs)

    def _init_id2lowerselect(self, objs):
        """Get all child item IDs for each item in dict keys."""
        for obj in objs:
            self._get_id2lowerselect(obj.item_id, obj)

    def _get_id2lowerselect(self, item_id, item_obj):
        """Add the child item IDs for one item object and their childs."""
        id2lowerselect = self.id2lowerselect
        if item_id in id2lowerselect:
            return id2lowerselect[item_id]
        child_ids = set()
        r2os = item_obj.relationship_rev
        sets = [r2os[r] for r in self.rset.intersection(r2os)]
        for child_obj in item_obj.children.union(*sets):
            child_id = child_obj.item_id
            child_ids.add(child_id)
            child_ids |= self._get_id2lowerselect(child_id, child_obj)
        id2lowerselect[item_id] = child_ids
        return child_ids

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
class CurNHigher:
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


# Copyright (C) 2010-present, DV Klopfenstein, H Tang, All rights reserved.
