"""Compare two or more sets of GO IDs. Best done using sections."""

from __future__ import print_function

__copyright__ = "Copyright (C) 2016-2019, DV Klopfenstein, H Tang. All rights reserved."
__author__ = "DV Klopfenstein"


from goatools.grouper.read_goids import read_sections
from goatools.grouper.grprdflts import GrouperDflts
from goatools.grouper.hdrgos import HdrgosSections
from goatools.grouper.grprobj import Grouper


# pylint: disable=too-few-public-methods
class Grouped(object):
    """Place all GO IDs in a Grouper."""

    kws_dict = set(['sections', 'S', 'go2nt', 'slims'])

    def __init__(self, gosubdag, godag_version, **kws):
        _kws = {k:v for k, v in kws.items() if k in self.kws_dict}
        self.grprdflt = GrouperDflts(gosubdag, _kws['slims'])
        self.ver_list = [godag_version, self.grprdflt.ver_goslims]
        self.sections = read_sections(self._get_secstr(**_kws), exclude_ungrouped=False)
        self.hdrobj = HdrgosSections(gosubdag, self.grprdflt.hdrgos_dflt, self.sections)
        _go2nt = _kws.get('go2nt')
        self.grprobj = Grouper("all", gosubdag.go_sources, self.hdrobj, gosubdag, go2nt=_go2nt)
        # for elem in self.grprobj.get_sections_2d():
        #     print(elem)
        # print('')

    @staticmethod
    def _get_secstr(**kws):
        """Return string containing sections file or sections module string."""
        if 'sections' in kws:
            return kws['sections']   # Sections text file
        if 'S' in kws:
            return kws['S']          # Sections module string


# Copyright (C) 2016-2019, DV Klopfenstein, H Tang. All rights reserved.
