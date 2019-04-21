"""Read annotations in GAF, GPAD, Entrez gene2go, or text format."""

__copyright__ = "Copyright (C) 2018-2019, DV Klopfenstein. All rights reserved."
__author__ = "DV Klopfenstein"

import sys
from goatools.anno.genetogo_reader import Gene2GoReader
from goatools.anno.gaf_reader import GafReader
from goatools.anno.gpad_reader import GpadReader
from goatools.anno.idtogos_reader import IdToGosReader


def get_objanno(fin_anno, anno_type=None, **kws):
    """Read annotations in GAF, GPAD, Entrez gene2go, or text format."""
    # kws get_objanno: taxids hdr_only prt allow_missing_symbol
    anno_type = get_anno_desc(fin_anno, anno_type)
    if anno_type is not None:
        if anno_type == 'gene2go':
            taxids = kws.get('taxids', None)
            return Gene2GoReader(fin_anno, taxids)
        if anno_type == 'gaf':
            hdr_only = kws.get('hdr_only', False)
            prt = kws.get('prt', sys.stdout)
            allow_missing_symbol = kws.get('allow_missing_symbol', False)
            return GafReader(fin_anno, hdr_only, prt, allow_missing_symbol)
        if anno_type == 'gpad':
            hdr_only = kws.get('hdr_only', False)
            return GpadReader(fin_anno, hdr_only)
        if anno_type == 'id2gos':
            return IdToGosReader(fin_anno)
    raise RuntimeError('UNEXPECTED ANNOTATION FILE FORMAT: {F} {D}'.format(
        F=fin_anno, D=anno_type))

def get_anno_desc(fin_anno, anno_type):
    """Indicate annotation format: gaf, gpad, NCBI gene2go, or id2gos."""
    if anno_type is not None:
        return anno_type
    if fin_anno[-7:] == 'gene2go':
        return 'gene2go'
    if fin_anno[-3:] == 'gaf':
        return 'gaf'
    if fin_anno[-3:] == 'gpa':
        return 'gpad'
    if fin_anno[-4:] == 'gpad':
        return 'gpad'


# Copyright (C) 2018-2019, DV Klopfenstein. All rights reserved
