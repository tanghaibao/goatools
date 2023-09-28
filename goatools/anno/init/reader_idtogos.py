"""Reads a Annotation File in text format with data in id2gos line"""

import timeit
import datetime
import collections as cx

from ...base import logger
from ...godag.consts import NAMESPACE2NS

__copyright__ = (
    "Copyright (C) 2016-present, DV Klopfenstein, H Tang. All rights reserved."
)
__author__ = "DV Klopfenstein"


# pylint: disable=too-few-public-methods
class InitAssc:
    """Initialize associations."""

    flds = ["DB_ID", "GO_ID"]

    def __init__(self, fin_anno, godag, namespaces, obsolete: str):
        tic = timeit.default_timer()
        self.godag = godag
        self.obsolete = obsolete
        self.id2gos = self._init_id2gos(fin_anno)
        self.nts = self.init_associations(namespaces)
        print(
            "HMS:{HMS} {N:7,} annotations READ: {ANNO} {NSs}".format(
                N=len(self.nts),
                ANNO=fin_anno,
                NSs=",".join(namespaces) if namespaces else "",
                HMS=str(datetime.timedelta(seconds=timeit.default_timer() - tic)),
            )
        )

    def init_associations(self, namespaces):
        """Get a list of namedtuples, one for each annotation."""
        nts = self._init_w_godag() if self.godag else self._init_dflt()
        if self.godag is None:
            if namespaces is not None:
                # pylint: disable=superfluous-parens
                logger.warning("GODAG NOT LOADED. IGNORING namespaces=%s", namespaces)
            return nts
        if namespaces == {"BP", "MF", "CC"}:
            return nts
        if not namespaces:
            return nts
        return [nt for nt in nts if nt.NS in namespaces]

    def _init_dflt(self):
        """Get a list of namedtuples, one for each annotation."""
        nts = []
        ntobj = cx.namedtuple("ntanno", self.flds)
        for itemid, gos in self.id2gos.items():
            for goid in gos:
                nts.append(ntobj(DB_ID=itemid, GO_ID=goid))
        return nts

    def _init_w_godag(self):
        """Get a list of namedtuples, one for each annotation."""
        nts = []
        ntobj = cx.namedtuple("ntanno", self.flds + ["NS"])
        s_godag = self.godag
        for itemid, gos in self.id2gos.items():
            to_add = set()
            for goid in gos:
                if goid not in s_godag:
                    logger.warning("%s NOT FOUND IN DAG", goid)
                    continue
                goobj = s_godag[goid]
                if goobj.is_obsolete:
                    if self.obsolete == "keep":
                        logger.warning("%s obsolete in DAG, kept", goid)
                        to_add.add(goid)
                    elif self.obsolete == "replace":
                        to_replace = set()
                        if "replaced_by" in goobj.__dict__ and goobj.replaced_by:
                            to_replace |= set(goobj.replaced_by.split(","))
                        if "consider" in goobj.__dict__ and goobj.consider:
                            to_replace |= goobj.consider
                        if to_replace:
                            logger.warning(
                                "%s obsolete in DAG, replaced by %s", goid, to_replace
                            )
                        else:
                            logger.warning("%s obsolete in DAG, no replacement", goid)
                        to_add |= to_replace
                    elif self.obsolete == "skip":
                        logger.warning("%s obsolete in DAG, skipped", goid)
                else:
                    to_add.add(goid)
            for goid in to_add:
                goobj = s_godag[goid]
                namespace = goobj.namespace
                nspc = NAMESPACE2NS.get(namespace, namespace) if goobj else ""
                nts.append(ntobj(DB_ID=itemid, GO_ID=goid, NS=nspc))
        return nts

    @staticmethod
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
        gene_int = None
        with open(assoc_fn) as ifstrm:
            for row in ifstrm:
                row = row.rstrip()
                atoms = row.split()
                if len(atoms) == 2:
                    gene_id, go_terms = atoms
                elif len(atoms) > 2 and row.count("\t") == 1:
                    gene_id, go_terms = row.split("\t")
                else:
                    continue
                gos = set(go_terms.split(";"))
                ## if no_top:
                ##     gos = gos.difference(top_terms)
                if gene_int is None:
                    gene_int = gene_id.isdigit()
                if gene_int:
                    gene_id = int(gene_id)
                assoc[gene_id] |= gos
        return assoc


# Copyright (C) 2016-present, DV Klopfenstein, H Tang. All rights reserved."
