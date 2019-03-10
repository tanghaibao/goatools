"""Manage GOATOOLS GOEA namedtuples."""

__copyright__ = "Copyright (C) 2010-2018, H Tang et al., All rights reserved."
__author__ = "DV Klopfenstein"

import collections as cx
from goatools.rpt.nts_xfrm import MgrNts


def get_study_items(goea_results):
    """Get all study items found in a GOATOOLS GOEA (e.g., geneids)."""
    return MgrNtGOEAs(goea_results).get_study_items()

def get_goea_nts_prt(goea_results, **kws):
    """Get namedtuples containing user-specified (or default) data from GOATOOLS GOEA results."""
    return MgrNtGOEAs(goea_results).get_goea_nts_prt(**kws)


class MgrNtGOEAs(object):
    """Manage GOATOOLS GOEA namedtuples."""

    def __init__(self, goea_results):
        self.goea_results = list(goea_results)

    def get_study_items(self):
        """Get all study items (e.g., geneids)."""
        study_items = set()
        for rec in self.goea_results:
            study_items |= rec.study_items
        return study_items

    def get_nts_strpval(self, fmt="{:8.2e}"):
        """Given GOEA namedtuples, return nts w/P-value in string format."""
        objntmgr = MgrNts(self.goea_results)
        dcts = objntmgr.init_dicts()
        # pylint: disable=line-too-long
        pval_flds = set(k for k in self._get_fieldnames(next(iter(self.goea_results))) if k[:2] == 'p_')
        for fld_float in pval_flds:
            fld_str = "s_" + fld_float[2:]
            objntmgr.add_f2str(dcts, fld_float, fld_str, fmt)
        return objntmgr.mknts(dcts)

    @staticmethod
    def dflt_sortby_objgoea(goea_res):
        """Default sorting of GOEA results."""
        return [getattr(goea_res, 'enrichment'),
                getattr(goea_res, 'namespace'),
                getattr(goea_res, 'p_uncorrected'),
                getattr(goea_res, 'depth'),
                getattr(goea_res, 'GO')]

    @staticmethod
    def dflt_sortby_ntgoea(ntgoea):
        """Default sorting of GOEA results stored in namedtuples."""
        return [ntgoea.enrichment,
                ntgoea.namespace,
                ntgoea.p_uncorrected,
                ntgoea.depth,
                ntgoea.GO]

    def get_goea_nts_prt(self, fldnames=None, **usr_kws):
        """Return list of namedtuples removing fields which are redundant or verbose."""
        kws = usr_kws.copy()
        if 'not_fldnames' not in kws:
            kws['not_fldnames'] = ['goterm', 'parents', 'children', 'id']
        if 'rpt_fmt' not in kws:
            kws['rpt_fmt'] = True
        return self.get_goea_nts_all(fldnames, **kws)

    def get_goea_nts_all(self, fldnames=None, **kws):
        """Get namedtuples containing user-specified (or default) data from GOEA results.

            Reformats data from GOEnrichmentRecord objects into lists of
            namedtuples so the generic table writers may be used.
        """
        # kws: prt_if indent itemid2name(study_items)
        data_nts = [] # A list of namedtuples containing GOEA results
        if not self.goea_results:
            return data_nts
        keep_if = kws.get('keep_if', None)
        rpt_fmt = kws.get('rpt_fmt', False)
        indent = kws.get('indent', False)
        # I. FIELD (column) NAMES
        not_fldnames = kws.get('not_fldnames', None)
        if fldnames is None:
            fldnames = self._get_fieldnames(self.goea_results[0])
        # Ia. Explicitly exclude specific fields from named tuple
        if not_fldnames is not None:
            fldnames = [f for f in fldnames if f not in not_fldnames]
        nttyp = cx.namedtuple("NtGoeaResults", " ".join(fldnames))
        goid_idx = fldnames.index("GO") if 'GO' in fldnames else None
        # II. Loop through GOEA results stored in a GOEnrichmentRecord object
        for goerec in self.goea_results:
            vals = self._get_field_values(goerec, fldnames, rpt_fmt, kws.get('itemid2name', None))
            if indent:
                vals[goid_idx] = "".join([goerec.get_indent_dots(), vals[goid_idx]])
            ntobj = nttyp._make(vals)
            if keep_if is None or keep_if(goerec):
                data_nts.append(ntobj)
        return data_nts

    @staticmethod
    def _get_field_values(item, fldnames, rpt_fmt=None, itemid2name=None):
        """Return fieldnames and values of either a namedtuple or GOEnrichmentRecord."""
        if hasattr(item, "_fldsdefprt"): # Is a GOEnrichmentRecord
            return item.get_field_values(fldnames, rpt_fmt, itemid2name)
        if hasattr(item, "_fields"): # Is a namedtuple
            return [getattr(item, f) for f in fldnames]

    @staticmethod
    def _get_fieldnames(item):
        """Return fieldnames of either a namedtuple or GOEnrichmentRecord."""
        if hasattr(item, "_fldsdefprt"): # Is a GOEnrichmentRecord
            return item.get_prtflds_all()
        if hasattr(item, "_fields"): # Is a namedtuple
            return item._fields

# Copyright (C) 2010-2018, H Tang et al., All rights reserved.
