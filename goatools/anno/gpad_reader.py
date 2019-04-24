"""Read a Gene Product Association Data (GPAD) and store the data in a Python object.
    Annotations available from the Gene Ontology Consortium:


    GPAD format:
        http://geneontology.org/page/gene-product-association-data-gpad-format
"""

import sys
import re
import collections as cx
from goatools.base import nopen
from goatools.anno.extensions.extensions import AnnotationExtensions
from goatools.anno.extensions.extension import AnnotationExtension
from goatools.anno.eco2group import ECO2GRP
from goatools.anno.annoreader_base import AnnoReaderBase
GET_DATE_YYYYMMDD = AnnoReaderBase.get_date_yyyymmdd

__copyright__ = "Copyright (C) 2016-2019, DV Klopfenstein, H Tang. All rights reserved."
__author__ = "DV Klopfenstein"


class GpadReader(AnnoReaderBase):
    """dRead a Gene Product Association Data (GPAD) and store the data in a Python object."""

    def __init__(self, filename=None, hdr_only=False):
        super(GpadReader, self).__init__('gpad', filename, hdr_only=hdr_only)
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
        ini = _InitAssc(fin_gpad)
        nts = ini.init_associations(hdr_only)
        self.hdr = ini.hdr
        return nts

# pylint: disable=too-few-public-methods
class _InitAssc(object):
    """Read annotation file and store a list of namedtuples."""

    # http://geneontology.org/page/gene-product-association-data-gpad-format
    gpadhdr = [ #              Col Req?     Cardinality    Example
        #                      --- -------- -------------- -----------------
        'DB',            #  0 required 1              UniProtKB
        'DB_ID',         #  1 required 1              P12345
        'Qualifier',     #  2 required 1 or greater   NOT
        'GO_ID',         #  3 required 1              GO:0003993
        # DB_Ref: set([''DOI:10.1002/sita.200600112', 'GO_REF:0000037', 'Reactome:R-HSA-6814682'])
        'DB_Reference',  #  4 required 1 or greater   set(['PMID:2676709',
        'ECO',           #  5 required 1              ECO:NNNNNNN
        'Evidence_Code', #    eco2group[ECO]
        'With_From',     #  6 optional 0 or greater   GO:0000346
        'Taxon',         #  7 optional 0 or 1         taxon:9606
        'Date',          #  8 required 1              20090118
        # Assigned_By: Ensembl FlyBase GO_Central GOC MGI Reactome UniProt WormBase
        'Assigned_By',   #  9 required 1              SGD
        # Annotations (Optional)
        'Extension',     # 10 optional 0 or greater
        'Properties',    # 11 optional 0 or greater
    ]

    gpad_columns = {"1.1" : gpadhdr}            # !gpad-version: 1.1

    # Expected numbers of columns for various versions
    exp_numcol = 12

    # Expected values for a Qualifier
    exp_qualifiers = set([
        'NOT', 'contributes_to', 'colocalizes_with', 'enables', 'involved_in',
        'part_of',
        # Seen starting 2019_03
        'is_active_in',
        # Seen starting 2018_09
        'acts_upstream_of',
        'acts_upstream_of_or_within',
        'acts_upstream_of_negative_effect',
        'acts_upstream_of_positive_effect',
        'acts_upstream_of_or_within_negative_effect',
        'acts_upstream_of_or_within_positive_effect',
    ])

    def __init__(self, filename):
        self.filename = filename
        self.hdr = None

    def _get_ntgpad(self, ntgpadobj, flds):
        """Convert fields from string to preferred format for GPAD ver 2.1 and 2.0."""
        is_set = False
        qualifiers = self._rd_fld_vals("Qualifier", flds[2], is_set)
        assert flds[3][:3] == 'GO:', 'UNRECOGNIZED GO({GO})'.format(GO=flds[3])
        db_reference = self._rd_fld_vals("DB_Reference", flds[4], is_set, 1)
        assert flds[5][:4] == 'ECO:', 'UNRECOGNIZED ECO({ECO})'.format(ECO=flds[3])
        with_from = self._rd_fld_vals("With_From", flds[6], is_set)
        taxons = self._get_taxon(flds[7])
        assert flds[8].isdigit(), 'UNRECOGNIZED DATE({D})'.format(D=flds[8])
        assert flds[9], '"Assigned By" VALUE WAS NOT FOUND'
        exten = self._get_extensions(flds[10])
        props = self._get_properties(flds[11])
        self._chk_qty_eq_1(flds, [0, 1, 3, 5, 8, 9])
        # Additional Formatting
        self._chk_qualifier(qualifiers)
        # Create list of values
        eco = flds[5]
        gpadvals = [
            flds[0],      #  0  DB
            flds[1],      #  1  DB_ID
            qualifiers,   #  3  Qualifier
            flds[3],      #  4  GO_ID
            db_reference, #  5  DB_Reference
            eco,          #  6  ECO
            ECO2GRP[eco],
            with_from,    #  7  With_From
            taxons,       # 12 Taxon
            GET_DATE_YYYYMMDD(flds[8]),      # 13 Date
            flds[9],      # 14 Assigned_By
            exten,        # 12 Extension
            props]        # 12 Annotation_Properties
        return ntgpadobj._make(gpadvals)

    @staticmethod
    def _rd_fld_vals(name, val, set_list_ft=True, qty_min=0, qty_max=None):
        """Further split a GPAD value within a single field."""
        if not val and qty_min == 0:
            return [] if set_list_ft else set()
        vals = val.split('|') # Use a pipe to separate entries
        num_vals = len(vals)
        assert num_vals >= qty_min, \
            "FLD({F}): MIN QUANTITY({Q}) WASN'T MET: {V}".format(
                F=name, Q=qty_min, V=vals)
        if qty_max is not None:
            assert num_vals <= qty_max, \
                "FLD({F}): MAX QUANTITY({Q}) EXCEEDED: {V}".format(
                    F=name, Q=qty_max, V=vals)
        return vals if set_list_ft else set(vals)

    @staticmethod
    def _get_taxon(taxon):
        """Return Interacting taxon ID | optional | 0 or 1 | gaf column 13."""
        if not taxon:
            return None
        ## assert taxon[:6] == 'taxon:', 'UNRECOGNIZED Taxon({Taxon})'.format(Taxon=taxon)
        ## taxid = taxon[6:]
        ## assert taxon[:10] == 'NCBITaxon:', 'UNRECOGNIZED Taxon({Taxon})'.format(Taxon=taxon)
        ## taxid = taxon[10:]
        # Get tzxon number: taxon:9606 NCBITaxon:9606
        sep = taxon.find(':')
        taxid = taxon[sep + 1:]
        assert taxid.isdigit(), "UNEXPECTED TAXON({T})".format(T=taxid)
        return int(taxid)

    def _get_properties(self, fldstr):
        """Return optional Annotation Properties (0 or greater)."""
        prop2val = {}
        props = self._rd_fld_vals("Properties", fldstr, False)  # Get set
        go_evidence = None
        for prop in props:
            # There can be more properties than 'go_evidence',
            # but currently we see only 'go_evidence'.
            # Upon encountering updates, evaluate and update code to support ...
            if prop[:12] == 'go_evidence=':
                assert go_evidence is None, "MORE THAN ONE EVIDENCE CODE FOUND"
                go_evidence = prop[12:]
            else:
                assert False, "UNPROGRAMMED PROPERTY({P})".format(P=prop)
        ## TBD: Is 'go_evidence' still used? Replaced by ECO? And eco2group
        ## assert go_evidence is not None, "go_evidence == None"
        ## prop2val['go_evidence'] = go_evidence
        if go_evidence is not None:
            prop2val['go_evidence'] = go_evidence
        return prop2val

    def _get_extensions(self, extline):
        """Return zero or greater Annotation Extensions, given a line of text."""
        # Extension examples:
        #   has_direct_input(UniProtKB:P37840),occurs_in(GO:0005576)
        #   part_of(UBERON:0006618),part_of(UBERON:0002302)
        #   occurs_in(CL:0000988)|occurs_in(CL:0001021)
        if not extline:
            return None
        exts = []
        for ext_lst in extline.split('|'):
            grp = []
            for ext in ext_lst.split(','):
                idx = ext.find('(')
                if idx != -1 and ext[-1] == ')':
                    grp.append(AnnotationExtension(ext[:idx], ext[idx+1:-1]))
                else:
                    # Ignore improperly formatted Extensions
                    sys.stdout.write('{F}: BAD Extension({E})\n'.format(F=self.filename, E=ext))
            exts.append(grp)
        return AnnotationExtensions(exts)

    # pylint: disable=too-many-locals
    def init_associations(self, hdr_only=False, prt=sys.stdout):
        """Read GPAD file. HTTP address okay. GZIPPED/BZIPPED file okay."""
        import timeit
        import datetime
        associations = []
        tic = timeit.default_timer()
        if self.filename is None:
            return associations
        ver = None
        ntgpadobj = None
        hdrobj = GpadHdr()
        ifstrm = nopen(self.filename)
        for lnum, line in enumerate(ifstrm, 1):
            # Read header
            if ntgpadobj is None:
                if line[0] == '!':
                    if ver is None and line[1:13] == 'gpa-version:':
                        ver = line[13:].strip()
                    hdrobj.chkaddhdr(line)
                else:
                    self.hdr = hdrobj.get_hdr()
                    if hdr_only:
                        return associations
                    ntgpadobj = cx.namedtuple("ntgpadobj", " ".join(self.gpad_columns[ver]))
            # Read data
            if ntgpadobj is not None:
                flds = self._split_line(line)
                try:
                    ntgpad = self._get_ntgpad(ntgpadobj, flds)
                    associations.append(ntgpad)
                except StandardError as inst:
                    import traceback
                    traceback.print_exc()
                    sys.stdout.write("\n  **FATAL: {MSG}\n\n".format(MSG=str(inst)))
                    sys.stdout.write("**FATAL: {FIN}[{LNUM}]:\n{L}\n".format(
                        FIN=self.filename, L=line, LNUM=lnum))
                    for idx, val in enumerate(flds):
                        sys.stdout.write('{I:2} {VAL}\n'.format(I=idx, VAL=val))
                    ## if datobj is not None:
                    ##     datobj.prt_line_detail(sys.stdout, line)
                    sys.exit(1)

        # GPAD file has been read
        prt.write('HMS:{HMS} {N:7,} annotations READ: {ANNO}\n'.format(
            N=len(associations), ANNO=self.filename,
            HMS=str(datetime.timedelta(seconds=(timeit.default_timer()-tic)))))
        return associations

    def _split_line(self, line):
        """Split line into field values."""
        line = line.rstrip('\r\n')
        flds = re.split('\t', line)
        assert len(flds) == self.exp_numcol, "EXPECTED({E}) COLUMNS, ACTUAL({A}): {L}".format(
            E=self.exp_numcol, A=len(flds), L=line)
        return flds

    def _chk_qualifier(self, qualifiers):
        """Check that qualifiers are expected values."""
        # http://geneontology.org/page/go-annotation-conventions#qual
        for qual in qualifiers:
            assert qual in self.exp_qualifiers, "UNEXPECTED QUALIFIER({Q}): {F}".format(
                Q=qual, F=self.filename)

    @staticmethod
    def _chk_qty_eq_1(flds, col_lst):
        """Check that these fields have only one value: required 1."""
        for col in col_lst:
            assert flds[col], "UNEXPECTED REQUIRED VALUE({V}) AT INDEX({R})".format(
                V=flds[col], R=col)


class GpadHdr(object):
    """Used to build a GPAD header."""

    cmpline = re.compile(r'^!(\w[\w\s-]+:.*)$')

    def __init__(self):
        self.gpadhdr = []

    def get_hdr(self):
        """Return GPAD header data as a string paragragh."""
        return "\n".join(self.gpadhdr)

    def chkaddhdr(self, line):
        """If this line contains desired header info, save it."""
        mtch = self.cmpline.search(line)
        if mtch:
            self.gpadhdr.append(mtch.group(1))

# Copyright (C) 2016-2019, DV Klopfenstein, H Tang. All rights reserved."
