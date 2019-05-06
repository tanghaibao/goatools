"""Read annotations in GAF, GPAD, Entrez gene2go, or text format."""

__copyright__ = "Copyright (C) 2018-2019, DV Klopfenstein. All rights reserved."
__author__ = "DV Klopfenstein"

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
            # kws: taxid taxids
            kws_gaf = {k:kws[k] for k in Gene2GoReader.exp_kws.intersection(kws.keys())}
            return Gene2GoReader(fin_anno, **kws)
        if anno_type == 'gaf':
            kws_gaf = {k:kws[k] for k in GafReader.exp_kws.intersection(kws.keys())}
            return GafReader(fin_anno, **kws_gaf)
        if anno_type == 'gpad':
            kws_gpad = {k:kws[k] for k in GpadReader.exp_kws.intersection(kws.keys())}
            return GpadReader(fin_anno, **kws_gpad)
        if anno_type == 'id2gos':
            kws_id2go = {k:kws[k] for k in IdToGosReader.exp_kws.intersection(kws.keys())}
            return IdToGosReader(fin_anno, **kws_id2go)
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
