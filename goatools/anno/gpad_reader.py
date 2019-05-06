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

    exp_kws = {'hdr_only', 'godag', 'namespaces'}

    def __init__(self, filename=None, **kws):
        super(GpadReader, self).__init__('gpad', filename,
                                         hdr_only=kws.get('hdr_only', False),
                                         godag=kws.get('godag'),
                                         namespaces=kws.get('namespaces'))
        self.qty = len(self.associations)

    def get_relation_cnt(self):
        """Return a Counter containing all relations contained in the Annotation Extensions."""
        ctr = cx.Counter()
        for ntgpad in self.associations:
            if ntgpad.Extension is not None:
                ctr += ntgpad.Extension.get_relations_cnt()
        return ctr

    def _init_associations(self, fin_gpad, **kws):
        """Read annotation file and store a list of namedtuples."""
        ini = InitAssc(fin_gpad, kws['godag'])
        nts = ini.init_associations(kws['hdr_only'], kws['namespaces'])
        self.hdr = ini.hdr
        return nts


# Copyright (C) 2016-2019, DV Klopfenstein, H Tang. All rights reserved."
