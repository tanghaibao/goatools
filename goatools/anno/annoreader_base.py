"""Read an Association File and store the data in a Python object."""

import sys
# import os
# import re
import collections as cx
# from datetime import datetime
from goatools.evidence_codes import EvidenceCodes

__copyright__ = "Copyright (C) 2016-2019, DV Klopfenstein, H Tang. All rights reserved."
__author__ = "DV Klopfenstein"


# pylint: disable=broad-except,too-few-public-methods,line-too-long
class AnnoReaderBase(object):
    """Reads a Gene Association File. Returns a Python object."""

    ## exp_kwdct = set(['allow_missing_symbol'])

    exp_qualifiers = set([
        # Seen in both GAF and gene2go
        'not', 'contributes_to', 'colocalizes_with',
        # Seen in gene2go
        # TBD: resolve gaf(not|contributes_to) v gene2go(contributes_to)
        'not contributes_to', 'not colocalizes_with',
    ])

    def __init__(self, filename=None, prt=sys.stdout):  ## , **kws):
        # kws: allow_missing_symbol
        ## self.kws = {k:v for k, v in kws.items() if k in self.exp_kwdct}
        self.filename = filename
        self.evobj = EvidenceCodes()
        self.associations = None

    ## def read_gaf(self, **kws):
    ##     """Read Gene Association File (GAF). Return data."""
    ##     # Simple associations
    ##     id2gos = cx.defaultdict(set)
    ##     # keyword arguments for choosing which GO IDs to keep
    ##     # Optional detailed associations split by taxid and having both ID2GOs & GO2IDs
    ##     taxid2asscs = kws.get('taxid2asscs', None)
    ##     b_geneid2gos = not kws.get('go2geneids', False)
    ##     evs = kws.get('evidence_set', None)
    ##     eval_nd = self._get_nd(kws.get('keep_ND', False))
    ##     eval_not = self._get_not(kws.get('keep_NOT', False))
    ##     # Optionally specify a subset of GOs based on their evidence.
    ##     # By default, return id2gos. User can cause go2geneids to be returned by:
    ##     #   >>> read_ncbi_gene2go(..., go2geneids=True
    ##     for ntgaf in self.associations:
    ##         if eval_nd(ntgaf) and eval_not(ntgaf):
    ##             if evs is None or ntgaf.Evidence_Code in evs:
    ##                 geneid = ntgaf.DB_ID
    ##                 go_id = ntgaf.GO_ID
    ##                 if b_geneid2gos:
    ##                     id2gos[geneid].add(go_id)
    ##                 else:
    ##                     id2gos[go_id].add(geneid)
    ##                 if taxid2asscs is not None:
    ##                     if ntgaf.Taxon:
    ##                         taxid = ntgaf.Taxon[0]
    ##                         taxid2asscs[taxid]['ID2GOs'][geneid].add(go_id)
    ##                         taxid2asscs[taxid]['GO2IDs'][go_id].add(geneid)
    ##     return id2gos # return simple associations

    ## @staticmethod
    ## def _get_nd(keep_nd):
    ##     """Allow GAF values always or never."""
    ##     if keep_nd:
    ##         return lambda nt: True
    ##     return lambda nt: nt.Evidence_Code != 'ND'

    ## @staticmethod
    ## def _get_not(keep_not):
    ##     """Allow GAF values always or never."""
    ##     if keep_not:
    ##         return lambda nt: True
    ##     return lambda nt: 'NOT' not in nt.Qualifier

    ## def _init_assn(self, fin_gaf, hdr_only, prt):
    ##     """Read GAF file. Store annotation data in a list of namedtuples."""
    ##     nts = self._read_gaf_nts(fin_gaf, hdr_only)
    ##     # GAF file has been read
    ##     if prt:
    ##         prt.write("  READ    {N:9,} associations: {FIN}\n".format(N=len(nts), FIN=fin_gaf))
    ##     # If there are illegal GAF lines ...
    ##     if self.datobj:
    ##         if self.datobj.ignored or self.datobj.illegal_lines:
    ##             self.datobj.prt_error_summary(fin_gaf)
    ##     return self.evobj.sort_nts(nts, 'Evidence_Code')

    ## def _read_gaf_nts(self, fin_gaf, hdr_only):
    ##     """Read GAF file. Store annotation data in a list of namedtuples."""
    ##     nts = []
    ##     ver = None
    ##     hdrobj = GafHdr()
    ##     datobj = None
    ##     lnum = line = -1
    ##     try:
    ##         with open(fin_gaf) as ifstrm:
    ##             for lnum, line in enumerate(ifstrm, 1):
    ##                 # Read header
    ##                 if datobj is None:
    ##                     if line[0] == '!':
    ##                         if ver is None and line[1:13] == 'gaf-version:':
    ##                             ver = line[13:].strip()
    ##                         hdrobj.chkaddhdr(line)
    ##                     else:
    ##                         self.hdr = hdrobj.get_hdr()
    ##                         if hdr_only:
    ##                             return nts
    ##                         datobj = GafData(ver, **self.kws)
    ##                 # Read data
    ##                 if datobj is not None and line[0] != '!':
    ##                     # print(lnum, line)
    ##                     ntgaf = datobj.get_ntgaf(line, lnum)
    ##                     if ntgaf is not None:
    ##                         nts.append(ntgaf)
    ##                     else:
    ##                         datobj.ignored.append((lnum, line))
    ##     except Exception as inst:
    ##         import traceback
    ##         traceback.print_exc()
    ##         sys.stderr.write("\n  **FATAL: {MSG}\n\n".format(MSG=str(inst)))
    ##         sys.stderr.write("**FATAL: {FIN}[{LNUM}]:\n{L}".format(FIN=fin_gaf, L=line, LNUM=lnum))
    ##         if datobj is not None:
    ##             datobj.prt_line_detail(sys.stdout, line)
    ##         sys.exit(1)
    ##     self.datobj = datobj
    ##     return nts

    def prt_summary_anno2ev(self, associations=None, prt=sys.stdout):
        """Print annotation/evidence code summary."""
        ctr = cx.Counter()
        if associations is None:
            associations = self.associations
        for ntanno in associations:
            evidence_code = ntanno.Evidence_Code
            if 'NOT' not in ntanno.Qualifier:
                ctr[evidence_code] += 1
            elif 'NOT' in ntanno.Qualifier:
                ctr["NOT {EV:3}".format(EV=ntanno.Evidence_Code)] += 1
            else:
                raise Exception("UNEXPECTED INFO")
        self.evobj.prt_ev_cnts(ctr, prt)



# Copyright (C) 2016-2019, DV Klopfenstein, H Tang. All rights reserved."
