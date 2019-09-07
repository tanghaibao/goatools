"""Reads a Annotation File in text format with data in id2gos line"""

import sys
from goatools.anno.annoreader_base import AnnoReaderBase
from goatools.anno.init.reader_idtogos import InitAssc

__copyright__ = "Copyright (C) 2016-2019, DV Klopfenstein, H Tang. All rights reserved."
__author__ = "DV Klopfenstein"


class IdToGosReader(AnnoReaderBase):
    """Reads a Annotation File in text format with data in id2gos line"""

    exp_kws = {'godag', 'namespaces'}

    def __init__(self, filename=None, **kws):
        self.id2gos = None  # ID to GO ID set as loaded from annotations file
        super(IdToGosReader, self).__init__('id2gos', filename,
                                            godag=kws.get('godag'),
                                            namespaces=kws.get('namespaces'))

    @staticmethod
    def wr_id2gos(fout_txt, id2gos):
        """Write annotations into a text file"""
        with open(fout_txt, 'w') as prt:
            for geneid, goset in sorted(id2gos.items()):
                prt.write('{GENE}\t{GOs}\n'.format(GENE=geneid, GOs=';'.join(sorted(goset))))
        print('  {N} annotations WROTE: {TXT}'.format(N=len(id2gos), TXT=fout_txt))

    def prt_summary_anno2ev(self, prt=sys.stdout):
        """Print a summary of all Evidence Codes seen in annotations"""
        prt.write('**NOTE: No evidence codes in associations: {F}\n'.format(F=self.filename))

    # pylint: disable=unused-argument
    def reduce_annotations(self, associations, options):
        """Return full annotations due to lack of Evidence_code or Qualifier in this format"""
        return associations

    def nts_ev_nd(self):
        """Get annotations where Evidence_code == 'ND' (No biological data)"""
        return []

    def nts_qual_not(self):
        """Get annotations having Qualifiers containing NOT"""
        return []

    def _init_associations(self, fin_anno, **kws):
        """Read annotation file and store a list of namedtuples."""
        ini = InitAssc(fin_anno, kws['godag'], kws['namespaces'])
        self.id2gos = ini.id2gos
        return ini.nts


# Copyright (C) 2016-2019, DV Klopfenstein, H Tang. All rights reserved."
