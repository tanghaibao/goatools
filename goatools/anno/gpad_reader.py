"""Read a Gene Product Association Data (GPAD) and store the data in a Python object.
    Annotations available from the Gene Ontology Consortium:


    GPAD format:
        http://geneontology.org/page/gene-product-association-data-gpad-format
"""

import collections as cx
from goatools.anno.annoreader_base import AnnoReaderBase
from goatools.anno.init.reader_gpad import InitAssc

__copyright__ = "Copyright (C) 2016-2019, DV Klopfenstein, H Tang. All rights reserved."
__author__ = "DV Klopfenstein"


class GpadReader(AnnoReaderBase):
    """dRead a Gene Product Association Data (GPAD) and store the data in a Python object."""

    def __init__(self, filename=None, hdr_only=False, godag=None):
        super(GpadReader, self).__init__('gpad', filename, godag, hdr_only=hdr_only)
        self.qty = len(self.associations)

    def get_relation_cnt(self):
        """Return a Counter containing all relations contained in the Annotation Extensions."""
        ctr = cx.Counter()
        for ntgpad in self.associations:
            if ntgpad.Extension is not None:
                ctr += ntgpad.Extension.get_relations_cnt()
        return ctr

    def _init_associations(self, fin_gpad, hdr_only=False):
        """Read annotation file and store a list of namedtuples."""
        ini = InitAssc(fin_gpad, self.godag)
        nts = ini.init_associations(hdr_only)
        self.hdr = ini.hdr
        return nts


# Copyright (C) 2016-2019, DV Klopfenstein, H Tang. All rights reserved."
