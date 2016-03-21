"""Read a GO Association File (GAF) and store in Python object.

    Annotations available from the Gene Ontology Consortium:
        http://geneontology.org/page/download-annotations

    GAF format:
        http://geneontology.org/page/go-annotation-file-formats
"""

import sys
import re
from collections import namedtuple
from base import nopen

__copyright__ = "Copyright (C) 2016, DV Klopfenstein, H Tang. All rights reserved."
__author__ = "DV Klopfenstein"

class GafReader(object):
    """Reads a Gene Annotation File (GAF). Returns a Python object."""

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

    gafhdr2 = [ #                   Col Req?     Cardinality    Example
        'Annotation_Extension', # 15 optional 0 or greater   part_of(CL:0000576)
        'Gene_Product_Form_ID', # 16 optional 0 or 1         UniProtKB:P12345-2
    ]

    gaf_columns = {
        "2.1" : gafhdr + gafhdr2, # !gaf-version: 2.1
        "2.0" : gafhdr + gafhdr2, # !gaf-version: 2.0
        "1.0" : gafhdr} # !gaf-version: 1.0

    # Expected numbers of columns for various versions
    gaf_numcol = {
        "2.1" : 17,
        "2.0" : 17,
        "1.0" : 15}

    # Expected values for a Qualifier
    exp_qualifiers = set(['NOT', 'contributes_to', 'colocalizes_with'])

    def __init__(self, filename=None, log=sys.stdout):
        self.filename = filename
        self.log = log
        self.associations = self.read_gaf(filename) if filename is not None else []

    def _get_ntgaf(self, ntgafobj, flds, ver):
        """Convert fields from string to preferred format for GAF ver 2.1 and 2.0."""
        # Cardinality
        is_set = False
        is_list = True
        qualifiers = self._rd_fld_vals("Qualifier", flds[3], is_set)
        db_reference = self._rd_fld_vals("DB_Reference", flds[5], is_set, 1)
        with_from = self._rd_fld_vals("With_From", flds[7], is_set)
        db_name = self._rd_fld_vals("DB_Name", flds[9], is_set, 0, 1)
        db_synonym = self._rd_fld_vals("DB_Synonym", flds[10], is_set)
        taxons = self._rd_fld_vals("Taxon", flds[12], is_list, 1, 2)
        self._chk_qty_eq_1(flds, [0, 1, 2, 4, 6, 8, 11, 13, 14])
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
        # Version 2.x had these additional fields not found in v1.0
        if ver[0] == '2':
            gafvals += [
                self._rd_fld_vals("Annotation_Extension", flds[15], is_set),
                self._rd_fld_vals("Gene_Product_Form_ID", flds[16], is_set)]
        return ntgafobj._make(gafvals)

    @staticmethod
    def _rd_fld_vals(name, val, set_list_ft=True, qty_min=0, qty_max=None):
        """Further split a GAF value within a single field."""
        if not val and qty_min == 0:
            return [] if set_list_ft else set()
        vals = val.split('|') # Use a pipe to separate entries
        num_vals = len(vals)
        assert num_vals >= qty_min, \
            "FLD({F}): MIN QUANTITY({Q}) NOT MET: {V}".format(
                F=name, Q=qty_min, V=vals)
        if qty_max is not None:
            assert num_vals <= qty_max, \
                "FLD({F}): MAX QUANTITY({Q}) EXCEEDED: {V}".format(
                    F=name, Q=qty_max, V=vals)
        return vals if set_list_ft else set(vals)

    def read_gaf(self, fin_gaf):
        """Read GAF file. HTTP address okay. GZIPPED/BZIPPED file okay."""
        ga_lst = []
        ifstrm = nopen(fin_gaf)
        ver = None
        ntgafobj = None
        exp_numcol = None
        for line in ifstrm:
            if ntgafobj is not None and not line.startswith('!'):
                flds = self._split_line(line, exp_numcol)
                ntgaf = self._get_ntgaf(ntgafobj, flds, ver)
                ga_lst.append(ntgaf)
            elif ntgafobj is None and line.startswith('!gaf-version:'):
                ver = line[13:].strip()
                ntgafobj = namedtuple("ntgafobj", " ".join(self.gaf_columns[ver]))
                exp_numcol = self.gaf_numcol[ver]
        self.log.write("  READ {N:,} items: {FIN}\n".format(N=len(ga_lst), FIN=fin_gaf))
        return ga_lst

    @staticmethod
    def _split_line(line, exp_numcol):
        """Split line into field values."""
        line = line.rstrip('\r\n')
        flds = re.split('\t', line)
        assert len(flds) == exp_numcol, "UNEXPECTED NUMBER OF COLUMNS"
        return flds

    def _chk_qualifier(self, qualifiers):
        """Check that qualifiers are expected values."""
        # http://geneontology.org/page/go-annotation-conventions#qual
        for qual in qualifiers:
            assert qual in self.exp_qualifiers, "UNEXPECTED QUALIFIER({Q})".format(Q=qual)

    @staticmethod
    def _chk_qty_eq_1(flds, col_lst):
        """Check that these fields have only one value: required 1."""
        for col in col_lst:
            assert flds[col], "UNEXPECTED REQUIRED VALUE({V}) AT INDEX({R})".format(
                V=flds[col], R=col)

    @staticmethod
    def _do_taxons(taxons):
        """Taxon"""
        taxons = [int(v[6:]) for v in taxons] # strip "taxon:"
        num_taxons = len(taxons)
        assert num_taxons == 1 or num_taxons == 2
        return taxons

# Copyright (C) 2016, DV Klopfenstein, H Tang. All rights reserved."
