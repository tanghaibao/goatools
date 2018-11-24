"""Sorts GO IDs or user-provided sections containing GO IDs."""

__copyright__ = "Copyright (C) 2016-2019, DV Klopfenstein, H Tang, All rights reserved."
__author__ = "DV Klopfenstein"

import collections as cx


class SorterGoIds(object):
    """Sorts grouped GO IDs.

       * Get a 2-D sorted list:
             goids = ["GO:HHHHHH0", [
                         "GO:UUUUU00",
                         ...
                         "GO:UUUUU0N"],
                     ["GO:HHHHHH1", [
                         "GO:UUUUU10",
                         ...
                         "GO:UUUUU1N"], ...

       * Get a flat sorted list of GO IDs
             goids = ["GO:HHHHHH0",
                      "GO:UUUUU00",
                      ...
                      "GO:UUUUU0N",
                      "GO:HHHHHH1",
                      "GO:UUUUU10",
                      ...
                      "GO:UUUUU1N" ...
    """

    def sortby(self, ntd):
        """Return function for sorting."""
        if 'reldepth' in self.grprobj.gosubdag.prt_attr['flds']:
            return [ntd.NS, -1*ntd.dcnt, ntd.reldepth]
        else:
            return [ntd.NS, -1*ntd.dcnt, ntd.depth]

    def __init__(self, grprobj, sortby=None, hdrgo_sortby=None):
        # User GO IDs grouped under header GO IDs are not sorted by the Grouper class.
        # Sort both user GO IDs in a group and header GO IDs across groups with these:
        #
        # H: hdrgo_sortby Sorts hdr GO IDs
        # U: sortby       Sorts user GO IDs
        # P: hdrgo_prt    If False, Removes GO IDs used as GO group headers; Leaves list in
        #                 sorted order, but removes header GO IDs which are not user GO IDs.
        #
        #  rm_h     hdr_sort      usr_sort       H  U  P
        #   ---     ------------  ------------   _  _  -
        #    NO     hdrgo_sortby  usrgo_sortby   H  U  T
        #   YES     hdrgo_sortby  usrgo_sortby   H  U  F
        self.grprobj = grprobj
        self.usrgo_sortby = self.sortby if sortby is None else sortby
        self.hdrgo_sortby = self._init_hdrgo_sortby(hdrgo_sortby, sortby)
        # Causes GO group headers to be removed in a flat list, if they are not user GO IDs.
        # Contains fields that can be used in sortby lambda keywords.
        assert grprobj is not None
        # print('SSSSSSSSSSS SorterGoIds(sortby={}, hdrgo_sortby={})'.format(sortby, hdrgo_sortby))
        # print('SSSSSSSSSSS SorterGoIds::self.usrgo_sortby:', self.usrgo_sortby)
        # print('SSSSSSSSSSS SorterGoIds::self.hdrgo_sortby:', self.hdrgo_sortby)

    # -- Sort header GO IDs and user GO IDs ----------------------------------------
    def get_nts_sorted(self, hdrgo_prt, hdrgos, hdrgo_sort):
        """Return a flat list of grouped and sorted GO terms."""
        nts_flat = []
        self.get_sorted_hdrgo2usrgos(hdrgos, nts_flat, hdrgo_prt, hdrgo_sort)
        return nts_flat

    def get_sorted_hdrgo2usrgos(self, hdrgos, flat_list=None, hdrgo_prt=True, hdrgo_sort=True):
        """Return GO IDs sorting using go2nt's namedtuple."""
        # Return user-specfied sort or default sort of header and user GO IDs
        sorted_hdrgos_usrgos = []
        h2u_get = self.grprobj.hdrgo2usrgos.get
        # Sort GO group headers using GO info in go2nt
        hdr_go2nt = self._get_go2nt(hdrgos)
        if hdrgo_sort is True:
            hdr_go2nt = sorted(hdr_go2nt.items(), key=lambda t: self.hdrgo_sortby(t[1]))
        for hdrgo_id, hdrgo_nt in hdr_go2nt:
            if flat_list is not None:
                if hdrgo_prt or hdrgo_id in self.grprobj.usrgos:
                    flat_list.append(hdrgo_nt)
            # Sort user GOs which are under the current GO header
            usrgos_unsorted = h2u_get(hdrgo_id)
            if usrgos_unsorted:
                usrgo2nt = self._get_go2nt(usrgos_unsorted)
                usrgont_sorted = sorted(usrgo2nt.items(), key=lambda t: self.usrgo_sortby(t[1]))
                usrgos_sorted, usrnts_sorted = zip(*usrgont_sorted)
                if flat_list is not None:
                    flat_list.extend(usrnts_sorted)
                sorted_hdrgos_usrgos.append((hdrgo_id, usrgos_sorted))
            else:
                sorted_hdrgos_usrgos.append((hdrgo_id, []))
        return cx.OrderedDict(sorted_hdrgos_usrgos)

    def _get_go2nt(self, goids):
        """Get go2nt for given goids."""
        go2nt_all = self.grprobj.go2nt
        return {go:go2nt_all[go] for go in goids}

    def _init_hdrgo_sortby(self, hdrgo_sortby, sortby):
        """Initialize header sort function."""
        if hdrgo_sortby is not None:
            return hdrgo_sortby
        if sortby is not None:
            return sortby
        return self.sortby


# Copyright (C) 2016-2019, DV Klopfenstein, H Tang, All rights reserved.
