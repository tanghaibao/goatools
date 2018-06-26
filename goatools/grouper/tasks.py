"""Tasks for grouping."""

__copyright__ = "Copyright (C) 2016-2018, DV Klopfenstein, H Tang, All rights reserved."
__author__ = "DV Klopfenstein"

from goatools.grouper.hdrgos import HdrgosSections


class SummarySec2dHdrGos(object):
    """Summary for a sections variable containing sets of header GO IDs."""

    def summarize_sec2hdrgos(self, sec2d_hdrgos):
        """Get counts of header GO IDs and sections."""
        hdrgos_all = set([])
        hdrgos_grouped = set()
        hdrgos_ungrouped = set()
        sections_grouped = set()
        for sectionname, hdrgos in sec2d_hdrgos:
            self._chk_hdrgoids(hdrgos)
            hdrgos_all.update(hdrgos)
            if sectionname != HdrgosSections.secdflt:
                hdrgos_grouped.update(hdrgos)
                sections_grouped.add(sectionname)
            else:
                hdrgos_ungrouped.update(hdrgos)
        return {'G': hdrgos_grouped,
                'S': sections_grouped,
                'U': hdrgos_all.difference(hdrgos_grouped)}

    def summarize_sec2hdrnts(self, sec2d_hdrnts):
        """Given namedtuples in each sectin, get counts of header GO IDs and sections."""
        sec2d_hdrgos = [(s, set(nt.GO for nt in nts)) for s, nts in sec2d_hdrnts]
        return self.summarize_sec2hdrgos(sec2d_hdrgos)

    @staticmethod
    def _chk_hdrgoids(hdrgos):
        """Check that hdrgo set is a set of GO IDs."""
        goid = next(iter(hdrgos))
        if isinstance(goid, str) and goid[:3] == "GO:":
            return
        assert False, "HDRGOS DO NOT CONTAIN GO IDs: {E}".format(E=goid)

# Copyright (C) 2016-2018, DV Klopfenstein, H Tang, All rights reserved.
