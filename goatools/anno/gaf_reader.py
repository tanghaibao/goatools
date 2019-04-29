"""Read a GO Annotation File (GAF) and store the data in a Python object.

    Annotations available from the Gene Ontology Consortium:
        http://geneontology.org/page/download-annotations

    GAF format:
        http://geneontology.org/page/go-annotation-file-formats
"""

import sys
from goatools.anno.annoreader_base import AnnoReaderBase
from goatools.anno.init.reader_gaf import GafData
from goatools.anno.init.reader_gaf import InitAssc

__copyright__ = "Copyright (C) 2016-2019, DV Klopfenstein, H Tang. All rights reserved."
__author__ = "DV Klopfenstein"


# pylint: disable=broad-except,too-few-public-methods,line-too-long
class GafReader(AnnoReaderBase):
    """Reads a Gene Annotation File (GAF). Returns a Python object."""

    def __init__(self, filename=None, hdr_only=False, **kws):
        super(GafReader, self).__init__(
            'gaf', filename, kws.get('godag'), hdr_only=hdr_only,
            prt=kws.get('prt', sys.stdout),
            allow_missing_symbol=kws.get('allow_missing_symbol', False))

    def read_gaf(self, **kws):
        """Read Gene Association File (GAF). Return associations."""
        return self.get_id2gos(**kws)

    def chk_associations(self, fout_err="gaf.err"):
        """Check that fields are legal in GAF"""
        obj = GafData("2.1")
        return obj.chk(self.associations, fout_err)

    def _init_associations(self, fin_gaf, hdr_only, prt, allow_missing_symbol):
        """Read annotation file and store a list of namedtuples."""
        ini = InitAssc()
        nts = ini.init_associations(fin_gaf, hdr_only, prt, allow_missing_symbol)
        self.hdr = ini.hdr
        return nts


# Copyright (C) 2016-2019, DV Klopfenstein, H Tang. All rights reserved."
