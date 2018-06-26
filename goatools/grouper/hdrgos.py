"""Manages GO group headers and optionally sections containing GO group headers."""

import collections as cx
from goatools.gosubdag.go_tasks import chk_goids

__copyright__ = "Copyright (C) 2016-2017, DV Klopfenstein, H Tang, All rights reserved."
__author__ = "DV Klopfenstein"


class HdrgosSections(object):
    """Manages GO group headers and optionally sections containing GO group headers."""

    secdflt = "Misc."

    def __init__(self, gosubdag, hdrgos_dflt, sections=None, **kws):
        self.gosubdag = gosubdag
        # sections is a 2-D list or is empty ([]) e.g, for section_name, hdrgos in sections:
        self.sections = self._init_sections(sections)
        self.hdrgos = self._init_hdrgos(hdrgos_dflt,
                                        kws.get('hdrgos', None), kws.get('add_dflt', True))
        self.hdrgo2sections = self._init_hdrgo2sections()
        assert len(set(self.hdrgo2sections.keys()).intersection(self.hdrgos)) == \
            len(self.hdrgo2sections)

    def get_sections(self, hdrgo, dflt_section=True):
        """Given a header GO, return the sections that contain it."""
        dflt_list = []
        # If the hdrgo is not in a section, return the default name for a section
        if dflt_section:
            dflt_list = [self.secdflt]
        return self.hdrgo2sections.get(hdrgo, dflt_list)

    def get_hdrgos(self):
        """Gets all possible GO grouping headers, regardardless if actually used or not."""
        return self.hdrgos.union(self.hdrgo2sections.keys())

    def num_hdrgos(self):
        """Number of GO group headers used for grouping."""
        return len(self.hdrgos) if self.hdrgos else 0

    def get_section_hdrgos(self):
        """Get the GO group headers explicitly listed in sections."""
        return set([h for _, hs in self.sections for h in hs]) if self.sections else set()

    @staticmethod
    def _chk_sections(sections):
        """Check format of user-provided 'sections' variable"""
        if sections:
            assert len(sections[0]) == 2, \
                "SECTIONS DATA MUST BE A 2-D LIST. FOUND: {S}".format(S=sections)
            for _, hdrgos in sections:
                chk_goids(hdrgos, "HdrgosSections::_chk_sections()")

    def _init_hdrgo2sections(self):
        """Return a dict with GO group headers as keys and section(s) as values."""
        hdrgo2sections = cx.defaultdict(list)
        for section_name, hdrgos in self.sections:
            for hdrgo in hdrgos:
                hdrgo2sections[hdrgo].append(section_name)
        return hdrgo2sections

    @staticmethod
    def _init_sections(sections):
        """Initialize sections."""
        if sections is None:
            return []
        return sections

    def _init_hdrgos(self, hdrgos_dflt, hdrgos_usr=None, add_dflt=True):
        """Initialize GO high"""
        # Use default GO group header values
        if (hdrgos_usr is None or hdrgos_usr is False) and not self.sections:
            return set(hdrgos_dflt)
        # Get GO group headers provided by user
        hdrgos_init = set()
        if hdrgos_usr:
            chk_goids(hdrgos_usr, "User-provided GO group headers")
            hdrgos_init |= set(hdrgos_usr)
        if self.sections:
            self._chk_sections(self.sections)
            hdrgos_sec = set([hg for _, hdrgos in self.sections for hg in hdrgos])
            chk_goids(hdrgos_sec, "User-provided GO group headers in sections")
            hdrgos_init |= hdrgos_sec
        # Add default depth-01 GOs to headers, if desired
        if add_dflt:
            return set(hdrgos_init).union(hdrgos_dflt)
        # Return user-provided GO grouping headers
        return hdrgos_init

# Copyright (C) 2016-2017, DV Klopfenstein, H Tang, All rights reserved.
