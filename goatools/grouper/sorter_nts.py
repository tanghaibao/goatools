"""Sorts GO IDs or user-provided sections containing GO IDs."""

__copyright__ = "Copyright (C) 2016-2019, DV Klopfenstein, H Tang, All rights reserved."
__author__ = "DV Klopfenstein"


class SorterNts(object):
    """Handles GO IDs in user-created sections.

       * Get a 2-D list of sections:
             sections = [
                 ['Immune', [
                      "GO:HHHHHH0", "GO:UUUUU00", ...  "GO:UUUUU0N", "GO:HHHHHH1", ...]],
                 ['Neuro', [
                      "GO:HHHHHH2", "GO:UUUUU20", ...  "GO:UUUUU2N", "GO:HHHHHH3", ...]],
             ]
       Also contains function for various tasks on grouped GO IDs:
         * Sort in various ways (sort by: p=value, depth, proximity to leaf-level, etc.):
           * Header GO ID groups
           * User GO IDs within a group
    """

    def __init__(self, sortgos, section_sortby=None):
        # User GO IDs grouped under header GO IDs are not sorted by the Grouper class.
        # Sort both user GO IDs in a group and header GO IDs across groups with these:
        # S: section_sortby (T=True, F=False, S=lambda sort function)
        # H: hdrgo_sortby Sorts hdr GO IDs
        # U: sortby       Sorts user GO IDs
        # P: hdrgo_prt    If True, Removes GO IDs used as GO group headers; Leaves list in
        #                 sorted order, but removes header GO IDs which are not user GO IDs.
        #
        #  rm_h     hdr_sort      usr_sort       S  H  U  P
        #   ---     ------------  ------------   _  _  _  -
        #    NO     hdrgo_sortby  usrgo_sortby   T  H  U  T
        #   YES     hdrgo_sortby  usrgo_sortby   T  H  U  F
        #    NO     section_order usrgo_sortby   F  -  U  T
        #   YES     section_order usrgo_sortby   F  -  U  F
        #   YES     |<----section_sortby---->|   S  -  -  -
        # print("SSSS SorterNts(sortgos, section_sortby={})".format(section_sortby))
        self.sortgos = sortgos  # SorterGoIds
        # section_sortby: True, False or None, or a sort_fnc
        self.section_sortby = section_sortby
        self.sections = self.sortgos.grprobj.hdrobj.sections
        # print('IIIIIIIIIIII SorterNts section_sortby', section_sortby)

    def get_sorted_nts_keep_section(self, hdrgo_prt):
        """Get 2-D list: 1st level is sections and 2nd level is grouped and sorted namedtuples."""
        section_nts = []
        # print("SSSS SorterNts:get_sorted_nts_keep_section(hdrgo_prt={})".format(hdrgo_prt))
        hdrgos_actual = self.sortgos.grprobj.get_hdrgos()
        hdrgos_secs = set()
        hdrgo_sort = False if self.section_sortby is False else True
        secname_dflt = self.sortgos.grprobj.hdrobj.secdflt
        for section_name, section_hdrgos_all in self.sections:
            #section_hdrgos_act = set(section_hdrgos_all).intersection(hdrgos_actual)
            section_hdrgos_act = [h for h in section_hdrgos_all if h in hdrgos_actual]
            hdrgos_secs |= set(section_hdrgos_act)
            nts_section = self.sortgos.get_nts_sorted(hdrgo_prt, section_hdrgos_act, hdrgo_sort)
            if nts_section:
                nts_section = self._get_sorted_section(nts_section)
                section_nts.append((section_name, nts_section))
        remaining_hdrgos = hdrgos_actual.difference(hdrgos_secs)
        # Add GO group headers not yet used under new section, Misc.
        if remaining_hdrgos:
            nts_section = self.sortgos.get_nts_sorted(hdrgo_prt, remaining_hdrgos, hdrgo_sort)
            if nts_section:
                nts_section = self._get_sorted_section(nts_section)
                section_nts.append((secname_dflt, nts_section))
        return section_nts

    def get_sorted_nts_omit_section(self, hdrgo_prt, hdrgo_sort):
        """Return a flat list of sections (wo/section names) with GO terms grouped and sorted."""
        nts_flat = []
        # print("SSSS SorterNts:get_sorted_nts_omit_section(hdrgo_prt={}, hdrgo_sort={})".format(
        #     hdrgo_prt, hdrgo_sort))
        hdrgos_seen = set()
        hdrgos_actual = self.sortgos.grprobj.get_hdrgos()
        for _, section_hdrgos_all in self.sections:
            #section_hdrgos_act = set(section_hdrgos_all).intersection(hdrgos_actual)
            section_hdrgos_act = [h for h in section_hdrgos_all if h in hdrgos_actual]
            hdrgos_seen |= set(section_hdrgos_act)
            self.sortgos.get_sorted_hdrgo2usrgos(
                section_hdrgos_act, nts_flat, hdrgo_prt, hdrgo_sort)
        remaining_hdrgos = set(self.sortgos.grprobj.get_hdrgos()).difference(hdrgos_seen)
        self.sortgos.get_sorted_hdrgo2usrgos(remaining_hdrgos, nts_flat, hdrgo_prt, hdrgo_sort)
        return nts_flat

    def _get_sorted_section(self, nts_section):
        """Sort GO IDs in each section, if requested by user."""
        #pylint: disable=unnecessary-lambda
        if self.section_sortby is True:
            return sorted(nts_section, key=lambda nt: self.sortgos.usrgo_sortby(nt))
        if self.section_sortby is False or self.section_sortby is None:
            return nts_section
        # print('SORT GO IDS IN A SECTION')
        return sorted(nts_section, key=lambda nt: self.section_sortby(nt))


# Copyright (C) 2016-2019, DV Klopfenstein, H Tang, All rights reserved.
