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
                        ntgaf = datobj.get_ntgaf(line)
                        if ntgaf is not None:
                            nts.append(ntgaf)
                        else:
                            self.prt_ignore_line(fin_gaf, line, lnum)
        except Exception as inst:
            import traceback
            traceback.print_exc()
            sys.stderr.write("\n  **FATAL in read_gaf: {MSG}\n\n".format(MSG=str(inst)))
            sys.stderr.write("**FATAL: {FIN}[{LNUM}]:\n{L}".format(FIN=fin_gaf, L=line, LNUM=lnum))
            sys.exit(1)
        # GAF file has been read
        if prt is not None:
            prt.write("  READ {N:,} associations: {FIN}\n".format(N=len(nts), FIN=fin_gaf))
        return self.evobj.sort_nts(nts, 'Evidence_Code')

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
    def prt_ignore_line(fin_gaf, line, lnum):
        """Print a message saying that we are ignoring an association line."""
        sys.stderr.write("**WARNING: IGNORED {FIN}[{LNUM}]:\n{L}\n\n".format(
            FIN=os.path.basename(fin_gaf), L=line, LNUM=lnum))

class GafData(object):
    """Extracts GAF fields from a GAF line."""

    spec_req1 = [0, 1, 2, 4, 6, 8, 11, 13, 14]

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
        self.exp_ncol = self.gaf_numcol[ver]
        self.req1 = self.spec_req1 if not allow_missing_symbol else [i for i in self.spec_req1 if i != 2]

    def get_ntgaf(self, line):
        """Return namedtuple filled with data."""
        flds = self._split_line(line)
        return self._get_ntgaf(flds)

    def _split_line(self, line):
        """Split line into field values."""
        line = line.rstrip('\r\n')
        flds = re.split('\t', line)
        assert len(flds) == self.exp_ncol, "EXP {E} COLUMNS, SAW {A}".format(
            A=len(flds), E=self.exp_ncol)
        return flds

    def _get_ntgaf(self, flds):
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
            gafvals += [
                self._rd_fld_vals("Annotation_Extension", flds[15], is_set),
                self._rd_fld_vals("Gene_Product_Form_ID", flds[16], is_set)]
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

    def _chk_qty_eq_1(self, flds):
        """Check that these fields have only one value: required 1."""
        for col in self.req1:
            if not flds[col]:
                sys.stderr.write("**ERROR: UNEXPECTED REQUIRED VAL({V}) FOR COL({R}):{H}: ".format(
                    V=flds[col], H=self.gafhdr[col], R=col))
                sys.stderr.write("{H0}({DB}) {H1}({ID})\n".format(
                    H0=self.gafhdr[0], DB=flds[0], H1=self.gafhdr[1], ID=flds[1]))
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
