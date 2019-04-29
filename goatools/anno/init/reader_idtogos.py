"""Reads a Annotation File in text format with data in id2gos line"""

import timeit
import datetime
import collections as cx
from goatools.godag.consts import Consts
NAMESPACE2NS = Consts.NAMESPACE2NS

__copyright__ = "Copyright (C) 2016-2019, DV Klopfenstein, H Tang. All rights reserved."
__author__ = "DV Klopfenstein"


# pylint: disable=too-few-public-methods
class InitAssc(object):
    """Initialize associations."""

    flds = ['DB_ID', 'GO_ID']

    def __init__(self, fin_anno, godag):
        tic = timeit.default_timer()
        self.godag = godag
        self.id2gos = self._init_id2gos(fin_anno)
        self.nts = self.init_associations()
        print('HMS:{HMS} {N:7,} annotations READ: {ANNO}'.format(
            N=len(self.nts), ANNO=fin_anno,
            HMS=str(datetime.timedelta(seconds=(timeit.default_timer()-tic)))))

    def init_associations(self):
        """Get a list of namedtuples, one for each annotation."""
        return self._init_w_godag() if self.godag else self._init_dflt()

    def _init_dflt(self):
        """Get a list of namedtuples, one for each annotation."""
        nts = []
        ntobj = cx.namedtuple('ntanno', self.flds)
        for itemid, gos in self.id2gos.items():
            for goid in gos:
                nts.append(ntobj(DB_ID=itemid, GO_ID=goid))
        return nts

    def _init_w_godag(self):
        """Get a list of namedtuples, one for each annotation."""
        nts = []
        ntobj = cx.namedtuple('ntanno', self.flds + ['NS'])
        for itemid, gos in self.id2gos.items():
            for goid in gos:
                goobj = self.godag.get(goid, '')
                nts.append(ntobj(
                    DB_ID=itemid,
                    GO_ID=goid,
                    NS=NAMESPACE2NS[goobj.namespace] if goobj else ''))
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
