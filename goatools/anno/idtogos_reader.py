"""Reads a Annotation File in text format with data in id2gos line"""

import sys
from goatools.anno.annoreader_base import AnnoReaderBase
from goatools.anno.init.reader_idtogos import InitAssc
from goatools.godag.consts import Consts
NAMESPACE2NS = Consts.NAMESPACE2NS

__copyright__ = "Copyright (C) 2016-2019, DV Klopfenstein, H Tang. All rights reserved."
__author__ = "DV Klopfenstein"


# pylint: disable=broad-except,too-few-public-methods,line-too-long
class IdToGosReader(AnnoReaderBase):
    """Reads a Annotation File in text format with data in id2gos line"""

    def __init__(self, filename=None, godag=None):
        self.id2gos = None
        super(IdToGosReader, self).__init__('id2gos', filename, godag)

    def prt_summary_anno2ev(self, prt=sys.stdout):
        """Print a summary of all Evidence Codes seen in annotations"""
        prt.write('**NOTE: No evidence codes in associations: {F}\n'.format(F=self.filename))

    def get_id2gos(self, **kws):
        """Return associations as a dict: id2gos"""
        return self._get_id2gos(self.associations, **kws) if kws else self.id2gos

    # pylint: disable=unused-argument
    def reduce_annotations(self, associations, options):
        """Return full annotations due to lack of Evidence_code or Qualifier in this format"""
        return self.associations

    def nts_ev_nd(self):
        """Get annotations where Evidence_code == 'ND' (No biological data)"""
        return []

    def nts_qual_not(self):
        """Get annotations having Qualifiers containing NOT"""
        return []

    def _init_associations(self, fin_anno):
        """Read annotation file and store a list of namedtuples."""
        ini = InitAssc(fin_anno, self.godag)
        # self.hdr = ini.flds
        self.id2gos = ini.id2gos
        return ini.nts


# Copyright (C) 2016-2019, DV Klopfenstein, H Tang. All rights reserved."
