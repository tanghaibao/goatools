"""Read a GO Association File (GAF) and store the data in a Python object.

    Annotations available from the Gene Ontology Consortium:
        http://geneontology.org/page/download-annotations

    GAF format:
        http://geneontology.org/page/go-annotation-file-formats
"""

import sys
import os
import re
import collections as cx
from goatools.evidence_codes import EvidenceCodes

__copyright__ = "Copyright (C) 2016-2018, DV Klopfenstein, H Tang. All rights reserved."
__author__ = "DV Klopfenstein"


# pylint: disable=broad-except,too-few-public-methods,line-too-long
class GafReader(object):
    """Reads a Gene Annotation File (GAF). Returns a Python object."""

    exp_kwdct = set(['allow_missing_symbol'])

    def __init__(self, filename=None, hdr_only=False, prt=sys.stdout, **kws):
        # kws: allow_missing_symbol
        self.kws = {k:v for k, v in kws.items() if k in self.exp_kwdct}
        self.filename = filename
        self.evobj = EvidenceCodes()
        # Initialize associations and header information
        self.hdr = None
        self.associations = self.read_gaf(filename, hdr_only, prt) if filename is not None else []

    def read_gaf(self, fin_gaf, hdr_only, prt):
        """Read GAF file. Store annotation data in a list of namedtuples."""
        nts = []
        ver = None
        hdrobj = GafHdr()
        datobj = None
        lnum = line = -1
        ignored = []
        try:
            with open(fin_gaf) as ifstrm:
                for lnum, line in enumerate(ifstrm, 1):
                    # Read header
                    if datobj is None:
                        if line[0] == '!':
                            if ver is None and line[1:13] == 'gaf-version:':
                                ver = line[13:].strip()
                            hdrobj.chkaddhdr(line)
                        else:
                            self.hdr = hdrobj.get_hdr()
                            if hdr_only:
                                return nts
                            datobj = GafData(ver, **self.kws)
                    # Read data
                    if datobj is not None and line[0] != '!':
                        # print(lnum, line)
                        ntgaf = datobj.get_ntgaf(line, lnum)
                        if ntgaf is not None:
                            nts.append(ntgaf)
                        else:
                            ignored.append((lnum, line))
        except Exception as inst:
            import traceback
            traceback.print_exc()
            sys.stderr.write("\n  **FATAL in read_gaf: {MSG}\n\n".format(MSG=str(inst)))
            sys.stderr.write("**FATAL: {FIN}[{LNUM}]:\n{L}".format(FIN=fin_gaf, L=line, LNUM=lnum))
            if datobj is not None:
                datobj.prt_line_detail(prt, line)
            sys.exit(1)
        # GAF file has been read
        self._prt_read_summary(prt, fin_gaf, nts, datobj, ignored)
        return self.evobj.sort_nts(nts, 'Evidence_Code')

    def _prt_read_summary(self, prt, fin_gaf, nts, datobj, ignored):
        """Print a summary about the GAF file that was read."""
        fout_log = self._prt_ignored_lines(ignored, datobj, fin_gaf) if ignored else None
        if prt is not None:
            prt.write("  READ    {N:9,} associations: {FIN}\n".format(N=len(nts), FIN=fin_gaf))
            if ignored:
                prt.write("  IGNORED {N:9,} associations: {FIN}\n".format(N=len(ignored), FIN=fout_log))

    def _prt_ignored_lines(self, ignored, datobj, fin_gaf):
        """Print ignored lines to a log file."""
        fout_log = "{}.log".format(fin_gaf)
        with open(fout_log, 'w') as prt:
            for lnum, line in ignored:
                self.prt_ignore_line(prt, fin_gaf, line, lnum)
                datobj.prt_line_detail(prt, line)
                prt.write("\n")
        return fout_log

    def prt_summary_anno2ev(self, prt=sys.stdout):
        """Print annotation/evidence code summary."""
        ctr = cx.Counter()
        for ntgaf in self.associations:
            evidence_code = ntgaf.Evidence_Code
            if 'NOT' not in ntgaf.Qualifier:
                ctr[evidence_code] += 1
            elif 'NOT' in ntgaf.Qualifier:
                ctr["NOT {EV}".format(EV=ntgaf.Evidence_Code)] += 1
            else:
                raise Exception("UNEXPECTED INFO")
        self.evobj.prt_ev_cnts(ctr, prt)

    @staticmethod
    def prt_ignore_line(prt, fin_gaf, line, lnum):
        """Print a message saying that we are ignoring an association line."""
        prt.write("**WARNING: BADLY FORMATTED LINE. IGNORED {FIN}[{LNUM}]:\n{L}\n".format(
            FIN=os.path.basename(fin_gaf), L=line, LNUM=lnum))

class GafData(object):
    """Extracts GAF fields from a GAF line."""

    spec_req1 = [0, 1, 2, 4, 6, 8, 11, 13, 14]

    req_str = ["REQ", "REQ", "REQ", "", "REQ", "REQ", "REQ", "", "REQ", "", "",
               "REQ", "REQ", "REQ", "REQ", "", ""]

    gafhdr = [ #           Col Req?     Cardinality    Example
        #                  --- -------- -------------- -----------------
        'DB',             #  0 required 1              UniProtKB
        'DB_ID',          #  1 required 1              P12345
        'DB_Symbol',      #  2 required 1              PHO3
        'Qualifier',      #  3 optional 0 or greater   NOT
        'GO_ID',          #  4 required 1              GO:0003993
        'DB_Reference',   #  5 required 1 or greater   PMID:2676709
        'Evidence_Code',  #  6 required 1              IMP
        'With_From',      #  7 optional 0 or greater   GO:0000346
        'Aspect',         #  8 required 1              F
        'DB_Name',        #  9 optional 0 or 1         Toll-like receptor 4
        'DB_Synonym',     # 10 optional 0 or greater   hToll|Tollbooth
        'DB_Type',        # 11 required 1              protein
        'Taxon',          # 12 required 1 or 2         taxon:9606
        'Date',           # 13 required 1              20090118
        'Assigned_By',    # 14 required 1              SGD
    ]

    #                            Col Required Cardinality  Example
    gafhdr2 = [ #                --- -------- ------------ -------------------
        'Annotation_Extension', # 15 optional 0 or greater part_of(CL:0000576)
        'Gene_Product_Form_ID', # 16 optional 0 or 1       UniProtKB:P12345-2
    ]

    gaf_columns = {
        "2.1" : gafhdr + gafhdr2, # !gaf-version: 2.1
        "2.0" : gafhdr + gafhdr2, # !gaf-version: 2.0
        "1.0" : gafhdr}           # !gaf-version: 1.0

    # Expected numbers of columns for various versions
    gaf_numcol = {
        "2.1" : 17,
        "2.0" : 17,
        "1.0" : 15}

    # Expected values for a Qualifier
    exp_qualifiers = set(['not', 'contributes_to', 'colocalizes_with'])

    def __init__(self, ver, allow_missing_symbol=False):
        self.ver = ver
        self.ntgafobj = cx.namedtuple("ntgafobj", " ".join(self.gaf_columns[ver]))
        self.req1 = self.spec_req1 if not allow_missing_symbol else [i for i in self.spec_req1 if i != 2]
        self.exp_mincol = 15  # Last required field is at the 15th column

    def get_ntgaf(self, line, lnum):
        """Return namedtuple filled with data."""
        flds = self.split_line(line)
        num_flds = len(flds)
        if num_flds >= self.exp_mincol:
            return self._get_ntgaf(flds, num_flds, lnum)

    @staticmethod
    def split_line(line):
        """Split line into field values."""
        line = re.split('\t', line)
        line[-1] = line[-1].rstrip('\r\n')
        return line

    def _get_ntgaf(self, flds, num_flds, lnum):
        """Convert fields from string to preferred format for GAF ver 2.1 and 2.0."""
        # Cardinality
        is_set = False
        is_list = True
        qualifiers = [t.lower() for t in self._rd_fld_vals("Qualifier", flds[3], is_set)]
        db_reference = self._rd_fld_vals("DB_Reference", flds[5], is_set, 1)
        with_from = self._rd_fld_vals("With_From", flds[7], is_set)
        db_name = self._rd_fld_vals("DB_Name", flds[9], is_set, 0)  # , 1)
        db_synonym = self._rd_fld_vals("DB_Synonym", flds[10], is_set)
        taxons = self._rd_fld_vals("Taxon", flds[12], is_list, 1, 2)
        if not self._chk_qty_eq_1(flds):
            return None
        # Additional Formatting
        taxons = self._do_taxons(taxons)
        self._chk_qualifier(qualifiers)
        # Create list of values
        gafvals = [
            flds[0],      # 0  DB
            flds[1],      # 1  DB_ID
            flds[2],      # 2  DB_Symbol
            qualifiers,   # 3  Qualifier
            flds[4],      # 4  GO_ID
            db_reference, # 5  DB_Reference
            flds[6],      # 6  Evidence_Code
            with_from,    # 7  With_From
            flds[8],      # 8  Aspect
            db_name,      # 9  DB_Name
            db_synonym,   # 10 DB_Synonym
            flds[11],     # 11 DB_Type
            taxons,       # 12 Taxon
            flds[12],     # 13 Date
            flds[13]]     # 14 Assigned_By
        # Version 2.x has these additional fields not found in v1.0
        if self.ver[0] == '2':
            # i=15) Annotation_Extension: optional 0 or greater; Ex: part_of(CL:0000576)
            if num_flds > 15:
                gafvals.append(self._rd_fld_vals("Annotation_Extension", flds[15], is_set))
            else:
                gafvals.append(None)
            # i=16) Gene_Product_Form_ID: optional 0 or 1;       Ex: UniProtKB:P12345-2
            if num_flds > 16:
                #self._prt_line_detail(sys.stdout, flds, lnum)
                gafvals.append(self._rd_fld_vals("Gene_Product_Form_ID", flds[16], is_set))
            else:
                gafvals.append(None)
        return self.ntgafobj._make(gafvals)

    @staticmethod
    def _rd_fld_vals(name, val, set_list_ft=True, qty_min=0, qty_max=None):
        """Further split a GAF value within a single field."""
        if not val and qty_min == 0:
            return [] if set_list_ft else set()
        vals = val.split('|') # Use a pipe to separate entries
        num_vals = len(vals)
        assert num_vals >= qty_min, \
            "FIELD({F}): MIN QUANTITY({Q}) WASN'T MET: {V}".format(F=name, Q=qty_min, V=vals)
        if qty_max is not None:
            assert num_vals <= qty_max, \
                "FIELD({F}): MAX QUANTITY({Q}) EXCEEDED: {V}".format(F=name, Q=qty_max, V=vals)
        return vals if set_list_ft else set(vals)

    def _chk_qualifier(self, qualifiers):
        """Check that qualifiers are expected values."""
        # http://geneontology.org/page/go-annotation-conventions#qual
        for qual in qualifiers:
            assert qual in self.exp_qualifiers, "UNEXPECTED QUALIFIER({Q})".format(Q=qual)

    def prt_line_detail(self, prt, line):
        """Print line header and values in a readable format."""
        values = self.split_line(line)
        self._prt_line_detail(prt, values)

    def _prt_line_detail(self, prt, values, lnum=""):
        """Print header and field values in a readable format."""
        data = zip(self.req_str, self.ntgafobj._fields, values)
        txt = ["{:2}) {:3} {:13} {}".format(i, req, hdr, val) for i, (req, hdr, val) in enumerate(data)]
        prt.write("{LNUM}\n{TXT}\n".format(LNUM=lnum, TXT="\n".join(txt)))

    def _chk_qty_eq_1(self, flds):
        """Check that these fields have only one value: required 1."""
        for col in self.req1:
            if not flds[col]:
                # sys.stderr.write("**ERROR: UNEXPECTED REQUIRED VAL({V}) FOR COL({R}):{H}: ".format(
                #     V=flds[col], H=self.gafhdr[col], R=col))
                # sys.stderr.write("{H0}({DB}) {H1}({ID})\n".format(
                #     H0=self.gafhdr[0], DB=flds[0], H1=self.gafhdr[1], ID=flds[1]))
                return False  # Check failed
        return True  # Check passed

    @staticmethod
    def _do_taxons(taxons):
        """Taxon"""
        taxons = [int(v[6:]) for v in taxons] # strip "taxon:"
        num_taxons = len(taxons)
        assert num_taxons == 1 or num_taxons == 2
        return taxons


class GafHdr(object):
    """Used to build a GAF header."""

    cmpline = re.compile(r'^!(\w[\w\s-]+:.*)$')

    def __init__(self):
        self.gafhdr = []

    def get_hdr(self):
        """Return GAF header data as a string paragragh."""
        return "\n".join(self.gafhdr)

    def chkaddhdr(self, line):
        """If this line contains desired header info, save it."""
        mtch = self.cmpline.search(line)
        if mtch:
            self.gafhdr.append(mtch.group(1))

# Copyright (C) 2016-2018, DV Klopfenstein, H Tang. All rights reserved."
