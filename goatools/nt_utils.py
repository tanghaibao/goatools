"""Utilities for combining namedtuples."""

__copyright__ = "Copyright (C) 2016-2018, DV Klopfenstein, H Tang. All rights reserved."
__author__ = "DV Klopfenstein"

import sys
import datetime
import collections as cx

def get_dict_w_id2nts(ids, id2nts, flds, dflt_null=""):
    """Return a new dict of namedtuples by combining "dicts" of namedtuples or objects."""
    assert len(ids) == len(set(ids)), "NOT ALL IDs ARE UNIQUE: {IDs}".format(IDs=ids)
    assert len(flds) == len(set(flds)), "DUPLICATE FIELDS: {IDs}".format(
        IDs=cx.Counter(flds).most_common())
    usr_id_nt = []
    # 1. Instantiate namedtuple object
    ntobj = cx.namedtuple("Nt", " ".join(flds))
    # 2. Fill dict with namedtuple objects for desired ids
    for item_id in ids:
        # 2a. Combine various namedtuples into a single namedtuple
        nts = [id2nt.get(item_id) for id2nt in id2nts]
        vals = _combine_nt_vals(nts, flds, dflt_null)
        usr_id_nt.append((item_id, ntobj._make(vals)))
    return cx.OrderedDict(usr_id_nt)

def get_list_w_id2nts(ids, id2nts, flds, dflt_null=""):
    """Return a new list of namedtuples by combining "dicts" of namedtuples or objects."""
    combined_nt_list = []
    # 1. Instantiate namedtuple object
    ntobj = cx.namedtuple("Nt", " ".join(flds))
    # 2. Fill dict with namedtuple objects for desired ids
    for item_id in ids:
        # 2a. Combine various namedtuples into a single namedtuple
        nts = [id2nt.get(item_id) for id2nt in id2nts]
        vals = _combine_nt_vals(nts, flds, dflt_null)
        combined_nt_list.append(ntobj._make(vals))
    return combined_nt_list

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

def wr_py_nts(fout_py, nts, docstring=None, varname="nts"):
    """Save namedtuples into a Python module."""
    if nts:
        with open(fout_py, 'w') as prt:
            prt.write('"""{DOCSTRING}"""\n\n'.format(DOCSTRING=docstring))
            prt.write("# Created: {DATE}\n".format(DATE=str(datetime.date.today())))
            prt_nts(prt, nts, varname)
            sys.stdout.write(" {N:7,} items WROTE: {PY}\n".format(N=len(nts), PY=fout_py))

def prt_nts(prt, nts, varname, spc='    '):
    """Print namedtuples into a Python module."""
    first_nt = nts[0]
    nt_name = type(first_nt).__name__
    prt.write("import collections as cx\n\n")
    prt.write("NT_FIELDS = [\n")
    for fld in first_nt._fields:
        prt.write('{SPC}"{F}",\n'.format(SPC=spc, F=fld))
    prt.write("]\n\n")
    prt.write('{NtName} = cx.namedtuple("{NtName}", " ".join(NT_FIELDS))\n\n'.format(
        NtName=nt_name))
    prt.write("# {N:,} items\n".format(N=len(nts)))
    prt.write("# pylint: disable=line-too-long\n")
    prt.write("{VARNAME} = [\n".format(VARNAME=varname))
    for ntup in nts:
        prt.write("{SPC}{NT},\n".format(SPC=spc, NT=ntup))
    prt.write("]\n")

def get_unique_fields(fld_lists):
    """Get unique namedtuple fields, despite potential duplicates in lists of fields."""
    flds = []
    fld_set = set([f for flst in fld_lists for f in flst])
    fld_seen = set()
    # Add unique fields to list of fields in order that they appear
    for fld_list in fld_lists:
        for fld in fld_list:
            # Add fields if the field has not yet been seen
            if fld not in fld_seen:
                flds.append(fld)
                fld_seen.add(fld)
    assert len(flds) == len(fld_set)
    return flds

# -- Internal methods ----------------------------------------------------------------
def _combine_nt_vals(lst0_lstn, flds, dflt_null):
    """Given a list of lists of nts, return a single namedtuple."""
    vals = []
    for fld in flds:
        fld_seen = False
        # Set field value using the **first** value seen in list of nt lists(lst0_lstn)
        for nt_curr in lst0_lstn:
            if hasattr(nt_curr, fld):
                vals.append(getattr(nt_curr, fld))
                fld_seen = True
                break
        # Set default value if GO ID or nt value is not present
        if fld_seen is False:
            vals.append(dflt_null)
    return vals

# Copyright (C) 2016-2018, DV Klopfenstein, H Tang. All rights reserved.
