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
        self.datobj = None
        self.associations = self._init_assn(filename, hdr_only, prt) if filename is not None else []

    def read_gaf(self, **kws):
        """Read Gene Association File (GAF). Return data."""
        # Simple associations
        id2gos = cx.defaultdict(set)
        # keyword arguments for choosing which GO IDs to keep
        # Optional detailed associations split by taxid and having both ID2GOs & GO2IDs
        taxid2asscs = kws.get('taxid2asscs', None)
        b_geneid2gos = not kws.get('go2geneids', False)
        evs = kws.get('evidence_set', None)
        eval_nd = self._get_nd(kws.get('keep_ND', False))
        eval_not = self._get_not(kws.get('keep_NOT', False))
        # Optionally specify a subset of GOs based on their evidence.
        # By default, return id2gos. User can cause go2geneids to be returned by:
        #   >>> read_ncbi_gene2go(..., go2geneids=True
        for ntgaf in self.associations:
            if eval_nd(ntgaf) and eval_not(ntgaf):
                if evs is None or ntgaf.Evidence_Code in evs:
                    geneid = ntgaf.DB_ID
                    go_id = ntgaf.GO_ID
                    if b_geneid2gos:
                        id2gos[geneid].add(go_id)
                    else:
                        id2gos[go_id].add(geneid)
                    if taxid2asscs is not None:
                        if ntgaf.Taxon:
                            taxid = ntgaf.Taxon[0]
                            taxid2asscs[taxid]['ID2GOs'][geneid].add(go_id)
                            taxid2asscs[taxid]['GO2IDs'][go_id].add(geneid)
        return id2gos # return simple associations

    @staticmethod
    def _get_nd(keep_nd):
        """Allow GAF values always or never."""
        if keep_nd:
            return lambda nt: True
        return lambda nt: nt.Evidence_Code != 'ND'

    @staticmethod
    def _get_not(keep_not):
        """Allow GAF values always or never."""
        if keep_not:
            return lambda nt: True
        return lambda nt: 'NOT' not in nt.Qualifier

    def _init_assn(self, fin_gaf, hdr_only, prt):
        """Read GAF file. Store annotation data in a list of namedtuples."""
        nts = self._read_gaf_nts(fin_gaf, hdr_only)
        # GAF file has been read
        if prt:
            prt.write("  READ    {N:9,} associations: {FIN}\n".format(N=len(nts), FIN=fin_gaf))
        # If there are illegal GAF lines ...
        if self.datobj:
            if self.datobj.ignored or self.datobj.illegal_lines:
                self.datobj.prt_error_summary(fin_gaf)
        return self.evobj.sort_nts(nts, 'Evidence_Code')

    def _read_gaf_nts(self, fin_gaf, hdr_only):
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
                        # print(lnum, line)
                        ntgaf = datobj.get_ntgaf(line, lnum)
                        if ntgaf is not None:
                            nts.append(ntgaf)
                        else:
                            datobj.ignored.append((lnum, line))
        except Exception as inst:
            import traceback
            traceback.print_exc()
            sys.stderr.write("\n  **FATAL: {MSG}\n\n".format(MSG=str(inst)))
            sys.stderr.write("**FATAL: {FIN}[{LNUM}]:\n{L}".format(FIN=fin_gaf, L=line, LNUM=lnum))
            if datobj is not None:
                datobj.prt_line_detail(sys.stdout, line)
            sys.exit(1)
        self.datobj = datobj
        return nts

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
    exp_qualifiers = set([
        'not', 'contributes_to', 'colocalizes_with',
        # Although thee go not appear in:
        #     http://geneontology.org/page/go-annotation-conventions#qual
        # they do appear in more than one July 2018 GAFs:
        #     'enables', 'involved_in', 'part_of', 
    ])

    def __init__(self, ver, allow_missing_symbol=False):
        self.ver = ver
        self.ntgafobj = cx.namedtuple("ntgafobj", " ".join(self.gaf_columns[ver]))
        self.req1 = self.spec_req1 if not allow_missing_symbol else [i for i in self.spec_req1 if i != 2]
        self.exp_mincol = 15  # Last required field is at the 15th column
        # Store information about illegal lines seen in a GAF file from the field
        self.ignored = []  # Illegal GAF lines that are ignored (e.g., missing an ID)
        self.illegal_lines = cx.defaultdict(list)  # GAF lines that are missing information (missing taxon)

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
        taxons = self._do_taxons(taxons, flds, lnum)
        self._chk_qualifier(qualifiers, flds, lnum)
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

    def _chk_qualifier(self, qualifiers, flds, lnum):
        """Check that qualifiers are expected values."""
        # http://geneontology.org/page/go-annotation-conventions#qual
        for qual in qualifiers:
            if qual not in self.exp_qualifiers:
                errname = 'UNEXPECTED QUALIFIER({QUAL})'.format(QUAL=qual)
                self.illegal_lines[errname].append((lnum, "\t".join(flds)))

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

    def _do_taxons(self, taxons, flds, lnum):
        """Taxon"""
        taxons_str = [v.split(':')[1] for v in taxons] # strip "taxon:"
        taxons_int = [int(s) for s in taxons_str if s]
        # taxons = [int(v[6:]) for v in taxons] # strip "taxon:"
        num_taxons = len(taxons_int)
        if taxons_int:
            assert num_taxons == 1 or num_taxons == 2
        else:
            self.illegal_lines['ILLEGAL TAXON'].append((lnum, "\t".join(flds)))
        return taxons_int

    def prt_error_summary(self, fin_gaf):
        """Print a summary about the GAF file that was read."""
        # Get summary of error types and their counts
        errcnts = []
        if self.ignored:
            errcnts.append("  {N:9,} IGNORED associations\n".format(N=len(self.ignored)))
        if self.illegal_lines:
            for err_name, errors in self.illegal_lines.items():
                errcnts.append("  {N:9,} {ERROR}\n".format(N=len(errors), ERROR=err_name))
        # Save error details into a log file
        fout_log = self._wrlog_details_illegal_gaf(fin_gaf, errcnts)
        sys.stdout.write("  WROTE GAF ERROR LOG: {LOG}:\n".format(LOG=fout_log))
        for err_cnt in errcnts:
            sys.stdout.write(err_cnt)

    def _wrlog_details_illegal_gaf(self, fin_gaf, err_cnts):
        """Print details regarding illegal GAF lines seen to a log file."""
        fout_log = "{}.log".format(fin_gaf)
        gaf_base = os.path.basename(fin_gaf)
        with open(fout_log, 'w') as prt:
            prt.write("ILLEGAL GAF ERROR SUMMARY:\n\n")
            for err_cnt in err_cnts:
                prt.write(err_cnt)
            prt.write("\n\nILLEGAL GAF ERROR DETAILS:\n\n")
            for lnum, line in self.ignored:
                prt.write("**WARNING: GAF LINE IGNORED: {FIN}[{LNUM}]:\n{L}\n".format(
                    FIN=gaf_base, L=line, LNUM=lnum))
                self.prt_line_detail(prt, line)
                prt.write("\n\n")
            for error, lines in self.illegal_lines.items():
                for lnum, line in lines:
                    prt.write("**WARNING: GAF LINE ILLEGAL({ERR}): {FIN}[{LNUM}]:\n{L}\n".format(
                        ERR=error, FIN=gaf_base, L=line, LNUM=lnum))
                    self.prt_line_detail(prt, line)
                    prt.write("\n\n")
        return fout_log


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
