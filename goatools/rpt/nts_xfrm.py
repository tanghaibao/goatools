"""Manage GOATOOLS GOEA namedtuples."""

__copyright__ = "Copyright (C) 2010-2018, H Tang et al., All rights reserved."
__author__ = "DV Klopfenstein"

import collections as cx


class MgrNts(object):
    """Manage GOATOOLS GOEA namedtuples."""

    def __init__(self, nts):
        self.nts = list(nts)

    def get_set(self, fieldname):
        """Get all study items (e.g., geneids)."""
        set_items = set()
        for ntdata in self.nts:
            set_items |= getattr(ntdata, fieldname)
        return set_items

    def mknts(self, add_dct):
        """Add information from add_dct to a new copy of namedtuples stored in nts."""
        nts = []
        assert len(add_dct) == len(self.nts)
        flds = list(next(iter(self.nts))._fields) + list(next(iter(add_dct)).keys())
        ntobj = cx.namedtuple("ntgoea", " ".join(flds))
        for dct_new, ntgoea in zip(add_dct, self.nts):
            dct_curr = ntgoea._asdict()
            for key, val in dct_new.items():
                dct_curr[key] = val
            nts.append(ntobj(**dct_curr))
        return nts

    def add_f2str(self, dcts, srcfld, dstfld, dstfmt):
        """Add a namedtuple field of type string generated from an existing namedtuple field."""
        # Example: f2str = objntmgr.add_f2str(dcts, "p_fdr_bh", "s_fdr_bh", "{:8.2e}")
        # ntobj = self.get_ntobj()
        # print(ntobj)
        assert len(dcts) == len(self.nts)
        for dct, ntgoea in zip(dcts, self.nts):
            valorig = getattr(ntgoea, srcfld)
            valstr = dstfmt.format(valorig)
            dct[dstfld] = valstr

    def get_ntobj(self):
        """Create namedtuple object with GOEA fields."""
        if self.nts:
            return cx.namedtuple("ntgoea", " ".join(vars(next(iter(self.nts))).keys()))

    def init_dicts(self):
        """Return a list of empty dicts to be filled with new data for revised namedtuples."""
        return [{} for _ in self.nts]


# Copyright (C) 2010-2018, H Tang et al., All rights reserved.
