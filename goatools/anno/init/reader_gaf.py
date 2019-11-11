"""Read a GO Annotation File (GAF) and store the data in a Python object.

    Annotations available from the Gene Ontology Consortium:
        http://current.geneontology.org/annotations/
"""

import sys
import os
import re
import collections as cx
import datetime
from goatools.anno.annoreader_base import AnnoReaderBase
from goatools.anno.init.utils import get_date_yyyymmdd
from goatools.anno.extensions.factory import get_extensions

__copyright__ = "Copyright (C) 2016-present, DV Klopfenstein, H Tang. All rights reserved."
__author__ = "DV Klopfenstein"


# pylint: disable=too-few-public-methods
class InitAssc:
    """Read annotation file and store a list of namedtuples."""

    def __init__(self, fin_gaf):
        self.fin_gaf = fin_gaf
        self.hdr = None
        self.datobj = None

    # pylint: disable=too-many-arguments
    def init_associations(self, hdr_only, prt, namespaces, allow_missing_symbol):
        """Read GAF file. Store annotation data in a list of namedtuples."""
        import timeit
        tic = timeit.default_timer()
        nts = self._read_gaf_nts(hdr_only, namespaces, allow_missing_symbol)
        # GAF file has been read
        if prt:
            prt.write('HMS:{HMS} {N:7,} annotations READ: {ANNO} {NSs}\n'.format(
                N=len(nts), ANNO=self.fin_gaf,
                NSs=','.join(namespaces) if namespaces else '',
                HMS=str(datetime.timedelta(seconds=(timeit.default_timer()-tic)))))
        # If there are illegal GAF lines ...
        if self.datobj:
            if self.datobj.ignored or self.datobj.illegal_lines:
                self.datobj.prt_error_summary(self.fin_gaf)
        return nts
        #### return self.evobj.sort_nts(nts, 'Evidence_Code')

    # pylint: disable=too-many-locals
    def _read_gaf_nts(self, hdr_only, namespaces, allow_missing_symbol):
        """Read GAF file. Store annotation data in a list of namedtuples."""
        nts = []
        ver = None
        hdrobj = GafHdr()
        datobj = None
        # pylint: disable=not-callable
        ntobj_make = None
        get_gafvals = None
        lnum = -1
        line = ''
        get_all_nss = namespaces is None or namespaces == {'BP', 'MF', 'CC'}
        try:
            with open(self.fin_gaf) as ifstrm:
                for lnum, line in enumerate(ifstrm, 1):
                    # Read data
                    if get_gafvals:
                        # print(lnum, line)
                        flds = line.split('\t')
                        nspc = GafData.aspect2ns[flds[8]]  # 8 GAF Aspect -> BP, MF, or CC
                        if get_all_nss or nspc in namespaces:
                            gafvals = get_gafvals(flds, nspc)
                            if gafvals:
                                nts.append(ntobj_make(gafvals))
                            else:
                                datobj.ignored.append((lnum, line))
                    # Read header
                    elif datobj is None:
                        if line[0] == '!':
                            if ver is None and line[1:13] == 'gaf-version:':
                                ver = line[13:].strip()
                            hdrobj.chkaddhdr(line)
                        else:
                            self.hdr = hdrobj.get_hdr()
                            if hdr_only:
                                return nts
                            datobj = GafData(ver, allow_missing_symbol)
                            get_gafvals = datobj.get_gafvals
                            ntobj_make = datobj.get_ntobj()._make
                            self._add_data0(nts, lnum, line, get_all_nss, namespaces, datobj)
        # pylint: disable=broad-except
        except Exception as inst:
            import traceback
            traceback.print_exc()
            sys.stderr.write("\n  **FATAL-gaf: {MSG}\n\n".format(MSG=str(inst)))
            sys.stderr.write("**FATAL-gaf: {FIN}[{LNUM}]:\n{L}".format(
                FIN=self.fin_gaf, L=line, LNUM=lnum))
            if datobj is not None:
                datobj.prt_line_detail(sys.stdout, line)
            sys.exit(1)
        self.datobj = datobj
        return nts

    @staticmethod
    def _add_data0(nts, lnum, line, get_all_nss, namespaces, datobj):
        """Do tasks upon finding the end of the header"""
        flds = line.split('\t')
        nspc = GafData.aspect2ns[flds[8]]  # 8 GAF Aspect -> BP, MF, or CC
        if get_all_nss or nspc in namespaces:
            gafvals = datobj.get_gafvals(flds, nspc)
            if gafvals:
                nts.append(datobj.get_ntobj()._make(gafvals))
            else:
                datobj.ignored.append((lnum, line))


class GafData:
    """Extracts GAF fields from a GAF line."""

    spec_req1 = [0, 1, 2, 4, 6, 8, 11, 13, 14]

    req_str = ["REQ", "REQ", "REQ", "", "REQ", "REQ", "REQ", "", "REQ", "", "",
               "REQ", "REQ", "REQ", "REQ", "", ""]

    aspect2ns = {'P':'BP', 'F':'MF', 'C':'CC'}

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
        'NS',             #  8 required 1              P->BP  F->MF  C->CC
        'DB_Name',        #  9 optional 0 or 1         Toll-like receptor 4
        'DB_Synonym',     # 10 optional 0 or greater   hToll|Tollbooth
        'DB_Type',        # 11 required 1              protein
        'Taxon',          # 12 required 1 or 2         taxon:9606
        'Date',           # 13 required 1              20090118
        'Assigned_By',    # 14 required 1              SGD
    ]

    #                            Col Required Cardinality  Example
    gafhdr2 = [ #                --- -------- ------------ -------------------
        'Extension', # 15 optional 0 or greater part_of(CL:0000576)
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

    def __init__(self, ver, allow_missing_symbol=False):
        self.ver = ver
        self.is_long = self._init_is_long(ver)
        self.flds = self.gaf_columns[self.ver]
        # pylint: disable=line-too-long
        self.req1 = self.spec_req1 if not allow_missing_symbol else [i for i in self.spec_req1 if i != 2]
        # Store information about illegal lines seen in a GAF file from the field
        self.ignored = []  # Illegal GAF lines that are ignored (e.g., missing an ID)
        self.illegal_lines = cx.defaultdict(list)  # GAF lines that are missing information (missing taxon)

    def _init_is_long(self, ver):
        """If the GAF version is 2.0 or 2.1, the GAF format is the long format (2 more cols)"""
        if ver is not None:
            return ver[0] == '2'
        print('\n**WARNING: NO VERSION LINE FOUND IN GAF FILE. USING:')
        print('!gaf-version: 2.1')
        self.ver = '2.1'
        return True

    @staticmethod
    def get_dfltdict():
        """Return default values for GAF data in a namedtuple"""
        return {
            'DB': 'database',        #  0 required 1              UniProtKB
            'DB_ID': 1,              #  1 required 1              P12345
            'DB_Symbol': 'Symbol1',  #  2 required 1              PHO3
            'Qualifier': set(),      #  3 optional 0 or greater   NOT
            'GO_ID': 'GO:0000001',   #  4 required 1              GO:0003993
            'DB_Reference': {'GO_REF:0000001',},   #  5 required 1 or greater   PMID:2676709
            'Evidence_Code': 'IDA',  #  6 required 1              IMP
            'With_From': set(),      #  7 optional 0 or greater   GO:0000346
            'NS': 'BP',              #  8 required 1              P->BP  F->MF  C->CC
            'DB_Name': set(),        #  9 optional 0 or 1         Toll-like receptor 4
            'DB_Synonym': set(),     # 10 optional 0 or greater   hToll|Tollbooth
            'DB_Type': 'protein',    # 11 required 1              protein
            'Taxon': [9606],         # 12 required 1 or 2         taxon:9606
            'Date': datetime.datetime.now().date(), # 13 required 1              20090118
            'Assigned_By': 'user',   # 14 required 1              SGD
            'Extension': set(),      # 15 optional 0 or greater part_of(CL:0000576)
            'Gene_Product_Form_ID':set(), # 16 optional 0 or 1       UniProtKB:P12345-2
        }

    def chk(self, annotations, fout_err):
        """Check annotations."""
        for idx, ntd in enumerate(annotations):
            self._chk_fld(ntd, "Qualifier")        # optional 0 or greater
            self._chk_fld(ntd, "DB_Reference", 1)  # required 1 or greater
            self._chk_fld(ntd, "With_From")        # optional 0 or greater
            self._chk_fld(ntd, "DB_Name", 0, 1)    # optional 0 or 1
            self._chk_fld(ntd, "DB_Synonym")       # optional 0 or greater
            self._chk_fld(ntd, "Taxon", 1, 2)
            flds = list(ntd)
            self._chk_qty_eq_1(flds)
            # self._chk_qualifier(ntd.Qualifier, flds, idx)
            if not ntd.Taxon or len(ntd.Taxon) not in {1, 2}:
                self.illegal_lines['BAD TAXON'].append(
                    (idx, '**{I}) TAXON: {NT}'.format(I=idx, NT=ntd)))
        if self.illegal_lines:
            self.prt_error_summary(fout_err)
        return not self.illegal_lines

    def get_ntobj(self):
        """Get namedtuple object specific to version"""
        return cx.namedtuple("ntgafobj", " ".join(self.flds))

    def get_gafvals(self, flds, nspc):
        """Convert fields from string to preferred format for GAF ver 2.1 and 2.0."""
        flds[3] = self._get_qualifier(flds[3])  # 3 Qualifier
        flds[5] = self._get_set(flds[5])     #  5 DB_Reference
        flds[7] = self._get_set(flds[7])     #  7 With_From
        flds[8] = nspc                       #  8 GAF Aspect field converted to BP, MF, or CC
        flds[9] = self._get_set(flds[9])     #  9 DB_Name
        flds[10] = self._get_set(flds[10])   # 10 DB_Synonym
        flds[12] = self._do_taxons(flds[12])   # 12 Taxon
        flds[13] = get_date_yyyymmdd(flds[13]) # 13 Date   20190406

        # Version 2.x has these additional fields not found in v1.0
        if self.is_long:
            flds[15] = get_extensions(flds[15])  # Extensions (or Annotation_Extension)
            flds[16] = self._get_set(flds[16].rstrip())
        else:
            flds[14] = self._get_set(flds[14].rstrip())
        return flds

    @staticmethod
    def _get_qualifier(val):
        """Get qualifiers. Correct for inconsistent capitalization in GAF files"""
        quals = set()
        if val == '':
            return quals
        for item in val.split('|'):
            item = item.lower()
            quals.add(item if item != 'not' else 'NOT')
        return quals

    @staticmethod
    def _get_set(val):
        """Further split a GAF value within a single field."""
        return set(val.split('|')) if val else set()

    @staticmethod
    def _get_list(val):
        """Further split a GAF value within a single field."""
        return val.split('|') if val else []

    def _chk_fld(self, ntd, name, qty_min=0, qty_max=None):
        """Further split a GAF value within a single field."""
        vals = getattr(ntd, name)
        num_vals = len(vals)
        if num_vals < qty_min:
            self.illegal_lines['MIN QTY'].append(
                (-1, "FIELD({F}): MIN QUANTITY({Q}) WASN'T MET: {V}".format(
                    F=name, Q=qty_min, V=vals)))
        if qty_max is not None:
            if num_vals > qty_max:
                pat = ("FIELD({F}): LINE({DB_ID} {GO_ID} {Evidence_Code}) "
                       "ERROR: MAX QUANTITY({Q}) EXCEEDED: {V}")
                err = pat.format(F=name, Q=qty_max, V=vals, **ntd._asdict())
                # TBD: DELETE THIS IF GAF FIXES THEIR FORMAT: -------------------------
                # https://github.com/geneontology/go-annotation/issues/2659
                # GAF files keep having errors. Try fixing and print warning if fixed
                # 2019_0921  1 ERROR IN: gramene_oryza.gaf
                vals.discard('')
                # 2019_0921 17 ERROR IN: pamgo_mgrisea.gaf  !Date Generated by GOC: 2019-07-01
                vals.discard('Similar to sp')
                vals.discard('Similar to tr')
                vals.discard('Similarities with tr')
                vals.discard('Similar to CA3529')
                vals.discard('Similar to CA3735')
                vals.discard('Similar to CA2345')
                # 2019_0921 5 ERROR IN: tair.gaf  !Date Generated by GOC: 2019-07-01
                vals.discard('Gb')
                # ----------------------------------------------------------------------
                if len(vals) > qty_max:
                    self.illegal_lines['MAX QTY'].append(
                        (-1, "{ERR}\n{NT}".format(ERR=err, NT=ntd)))
                else:
                    print('**WARNING GAF FILE: {ERR}'.format(ERR=err))

    def _chk_qualifier(self, qualifiers, flds, lnum):
        """Check that qualifiers are expected values."""
        # http://geneontology.org/page/go-annotation-conventions#qual
        for qual in qualifiers:
            if qual not in AnnoReaderBase.exp_qualifiers:
                errname = 'UNEXPECTED QUALIFIER({QUAL})'.format(QUAL=qual)
                self.illegal_lines[errname].append((lnum, "\t".join(flds)))

    def prt_line_detail(self, prt, line):
        """Print line header and values in a readable format."""
        values = line.split('\t')
        self._prt_line_detail(prt, values)

    def _prt_line_detail(self, prt, values, lnum=""):
        """Print header and field values in a readable format."""
        data = zip(self.req_str, self.flds, values)
        pat = "{:2}) {:3} {:20} {}"
        txt = [pat.format(i, req, hdr, val) for i, (req, hdr, val) in enumerate(data)]
        prt.write("{LNUM}\n{TXT}\n".format(LNUM=lnum, TXT="\n".join(txt)))

    def _chk_qty_eq_1(self, flds):
        """Check that these fields have only one value: required 1."""
        for col in self.req1:
            if not flds[col]:
                self.illegal_lines['QTY 1'].append(
                    (-1, "**ERROR: UNEXPECTED REQUIRED VAL({V}) FOR COL({R}):{H}: ".format(
                        V=flds[col], H=self.gafhdr[col], R=col)))
                self.illegal_lines['QTY 1'].append((-1, "{H0}({DB}) {H1}({ID})\n".format(
                    H0=self.gafhdr[0], DB=flds[0], H1=self.gafhdr[1], ID=flds[1])))

    def _do_taxons(self, taxon_str):
        """Taxon"""
        taxons = self._get_list(taxon_str)
        taxons_str = [v.split(':')[1] for v in taxons] # strip "taxon:"
        taxons_int = [int(s) for s in taxons_str if s]
        return taxons_int

    def prt_error_summary(self, fout_err):
        """Print a summary about the GAF file that was read."""
        # Get summary of error types and their counts
        errcnts = []
        if self.ignored:
            errcnts.append("  {N:9,} IGNORED associations\n".format(N=len(self.ignored)))
        if self.illegal_lines:
            for err_name, errors in self.illegal_lines.items():
                errcnts.append("  {N:9,} {ERROR}\n".format(N=len(errors), ERROR=err_name))
        # Save error details into a log file
        fout_log = self._wrlog_details_illegal_gaf(fout_err, errcnts)
        sys.stdout.write("  WROTE GAF ERROR LOG: {LOG}:\n".format(LOG=fout_log))
        for err_cnt in errcnts:
            sys.stdout.write(err_cnt)

    def _wrlog_details_illegal_gaf(self, fout_err, err_cnts):
        """Print details regarding illegal GAF lines seen to a log file."""
        # fout_err = "{}.log".format(fin_gaf)
        gaf_base = os.path.basename(fout_err)
        with open(fout_err, 'w') as prt:
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
        return fout_err


class GafHdr:
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


# Copyright (C) 2016-present, DV Klopfenstein, H Tang. All rights reserved."
