"""Reads a Annotation File in text format with data in id2gos line"""

import sys
import collections as cx
from goatools.anno.annoreader_base import AnnoReaderBase

__copyright__ = "Copyright (C) 2016-2019, DV Klopfenstein, H Tang. All rights reserved."
__author__ = "DV Klopfenstein"


# pylint: disable=broad-except,too-few-public-methods,line-too-long
class IdToGosReader(AnnoReaderBase):
    """Reads a Annotation File in text format with data in id2gos line"""

    def __init__(self, filename=None):  # , **kws):
        self.id2gos = None
        super(IdToGosReader, self).__init__('id2gos', filename)

    def prt_summary_anno2ev(self, prt=sys.stdout):
        """Print a summary of all Evidence Codes seen in annotations"""
        prt.write('**NOTE: No evidence codes in associations: {F}\n'.format(F=self.filename))

    def get_id2gos(self, **kws):
        """Return associations as a dict: id2gos"""
        return self.id2gos

    # pylint: disable=unused-argument
    def reduce_annotations(self, associations, options):
        """Return full annotations due to lack of Evidence_code or Qualifier in this format"""
        return self.associations

    # - initialization -------------------------------------------------------------------------
    def _init_associations(self, fin_anno):
        """Read annotation file and store a list of namedtuples."""
        ini = _InitAssc(fin_anno)
        # self.hdr = ini.flds
        self.id2gos = ini.id2gos
        return ini.nts

class _InitAssc(object):

    flds = ['DB_ID', 'GO_ID']

    def __init__(self, fin_anno):
        import timeit
        import datetime
        tic = timeit.default_timer()
        self.id2gos = self._init_id2gos(fin_anno)
        self.nts = self.init_associations()
        print('HMS:{HMS} {N:7,} annotations READ: {ANNO}'.format(
            N=len(self.nts), ANNO=fin_anno,
            HMS=str(datetime.timedelta(seconds=(timeit.default_timer()-tic)))))

    def init_associations(self):
        """Get a list of namedtuples, one for each annotation."""
        nts = []
        ntobj = cx.namedtuple('ntanno', self.flds)
        for itemid, gos in self.id2gos.items():
            for goid in gos:
                nts.append(ntobj(DB_ID=itemid, GO_ID=goid))
        return nts

    @staticmethod
    #### def read_associations(assoc_fn, no_top=False):
    def _init_id2gos(assoc_fn):  ##, no_top=False):
        """
        Reads a gene id go term association file. The format of the file
        is as follows:

        AAR1	GO:0005575;GO:0003674;GO:0006970;GO:0006970;GO:0040029
        AAR2	GO:0005575;GO:0003674;GO:0040029;GO:0009845
        ACD5	GO:0005575;GO:0003674;GO:0008219
        ACL1	GO:0005575;GO:0003674;GO:0009965;GO:0010073
        ACL2	GO:0005575;GO:0003674;GO:0009826
        ACL3	GO:0005575;GO:0003674;GO:0009826;GO:0009965

        Also, the following format is accepted (gene ids are repeated):

        AAR1	GO:0005575
        AAR1    GO:0003674
        AAR1    GO:0006970
        AAR2	GO:0005575
        AAR2    GO:0003674
        AAR2    GO:0040029

        :param assoc_fn: file name of the association
        :return: dictionary having keys: gene id, values set of GO terms
        """
        assoc = cx.defaultdict(set)
        ## top_terms = set(['GO:0008150', 'GO:0003674', 'GO:0005575']) # BP, MF, CC
        for row in open(assoc_fn, 'r'):
            atoms = row.split()
            if len(atoms) == 2:
                gene_id, go_terms = atoms
            elif len(atoms) > 2 and row.count('\t') == 1:
                gene_id, go_terms = row.split("\t")
            else:
                continue
            gos = set(go_terms.split(";"))
            ## if no_top:
            ##     gos = gos.difference(top_terms)
            assoc[gene_id] |= gos
        return assoc


# Copyright (C) 2016-2019, DV Klopfenstein, H Tang. All rights reserved."
