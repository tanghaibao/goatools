"""Reads a Annotation File in text format with data in id2gos line"""

import sys

from ..base import logger
from .annoreader_base import AnnoReaderBase
from .init.reader_idtogos import InitAssc, _parse_go_terms

__copyright__ = "Copyright (C) 2016-2019, DV Klopfenstein, H Tang. All rights reserved."
__author__ = "DV Klopfenstein"


class IdToGosReader(AnnoReaderBase):
    """Reads a Annotation File in text format with data in id2gos line"""

    exp_kws = {"godag", "namespaces", "obsolete"}

    def __init__(self, filename=None, **kws):
        self.id2gos = None  # ID to GO ID set as loaded from annotations file
        super().__init__(
            "id2gos",
            filename,
            godag=kws.get("godag"),
            namespaces=kws.get("namespaces"),
            obsolete=kws.get("obsolete"),
        )

    @staticmethod
    def wr_id2gos(fout_txt, id2gos):
        """Write annotations into a text file"""
        with open(fout_txt, "w") as prt:
            for geneid, goset in sorted(id2gos.items()):
                prt.write(
                    "{GENE}\t{GOs}\n".format(GENE=geneid, GOs=";".join(sorted(goset)))
                )
        print(f"  {len(id2gos)} annotations WROTE: {fout_txt}")

    def prt_summary_anno2ev(self, prt=sys.stdout):
        """Print a summary of all Evidence Codes seen in annotations"""
        logger.info("No evidence codes in associations: %s", self.filename)

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

    def chk_associations(self, fout_err=None):
        """Check GO IDs in the id2gos association file for validity.

        Scans each GO ID in the file and reports any that do not match the
        expected format GO:NNNNNNN (7 digits).  Call this after loading
        annotations to identify formatting problems such as extra spaces
        introduced by using "; " instead of ";" as a separator.

        :param fout_err: optional filename to write invalid-ID lines to
        :return: number of invalid GO IDs found
        """
        num_invalid = 0
        prt_err = open(fout_err, "w") if fout_err else None
        try:
            with open(self.filename) as ifstrm:
                for lineno, row in enumerate(ifstrm, 1):
                    row = row.rstrip()
                    atoms = row.split()
                    if len(atoms) == 2:
                        _, go_terms = atoms
                    elif len(atoms) > 2 and row.count("\t") == 1:
                        _, go_terms = row.split("\t")
                    else:
                        continue
                    _, invalid = _parse_go_terms(go_terms)
                    for bad in invalid:
                        num_invalid += 1
                        logger.warning(
                            "LINE %d: INVALID GO ID '%s' IN: %s",
                            lineno,
                            bad,
                            self.filename,
                        )
                        if prt_err:
                            prt_err.write(
                                "LINE {L}: INVALID GO ID '{G}'\n".format(
                                    L=lineno, G=bad
                                )
                            )
        finally:
            if prt_err:
                prt_err.close()
        return num_invalid

    def _init_associations(self, fin_anno, **kws):
        """Read annotation file and store a list of namedtuples."""
        ini = InitAssc(fin_anno, kws["godag"], kws["namespaces"], kws["obsolete"])
        self.id2gos = ini.id2gos
        return ini.nts


# Copyright (C) 2016-2019, DV Klopfenstein, H Tang. All rights reserved."
