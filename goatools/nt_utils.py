"""Utilities for combining namedtuples."""

import collections as cx

# -- combine_id2nts ------------------------------------------------------------------
def combine_id2nts(ids, id2nts, flds, dflt_null=""):
    """Return a new dict of namedtuples by combining "dicts" of namedtuples or objects."""
    combined_nt_go2nt = []
    # 1. Instantiate namedtuple object
    ntobj = cx.namedtuple("Nt", " ".join(flds))
    # 2. Fill dict with namedtuple objects for desired ids
    for item_id in ids:
        # 2a. Combine various namedtuples into a single namedtuple
        nts = [id2nt.get(item_id) for id2nt in id2nts]
        vals = _combine_nt_vals(nts, flds, dflt_null)
        combined_nt_go2nt.append(ntobj._make(vals))
    return combined_nt_go2nt

# -- combine_nt_lists ----------------------------------------------------------------
def combine_nt_lists(lists, flds, dflt_null=""):
    """Return a new list of namedtuples by zipping "lists" of namedtuples or objects."""
    combined_nt_list = []
    # Check that all lists are the same length
    lens = [len(lst) for lst in lists]
    assert len(set(lens)) == 1, \
        "LIST LENGTHS MUST BE EQUAL: {Ls}".format(Ls=" ".join(str(l) for l in lens))
    # 1. Instantiate namedtuple object
    ntobj = cx.namedtuple("Nt", " ".join(flds))
    # 2. Loop through zipped list
    for lst0_lstn in zip(*lists):
        # 2a. Combine various namedtuples into a single namedtuple
        combined_nt_list.append(ntobj._make(_combine_nt_vals(lst0_lstn, flds, dflt_null)))
    return combined_nt_list

def _combine_nt_vals(lst0_lstn, flds, dflt_null):
    """Given a list of lists of nts, return a single namedtuple."""
    vals = []
    for fld in flds:
        fld_seen = False
        for nt_curr in lst0_lstn:
            if hasattr(nt_curr, fld):
                vals.append(getattr(nt_curr, fld))
                fld_seen = True
                break
        if fld_seen is False:
            vals.append(dflt_null) # Default val if GO id or nt val is not present
    return vals
