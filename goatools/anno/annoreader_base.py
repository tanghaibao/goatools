"""Read an Association File and store the data in a Python object."""

import sys
import timeit
import datetime
from datetime import date
# import os
# import re
import collections as cx
from goatools.evidence_codes import EvidenceCodes
from goatools.anno.opts import AnnoOptions

__copyright__ = "Copyright (C) 2016-2019, DV Klopfenstein, H Tang. All rights reserved."
__author__ = "DV Klopfenstein"


# pylint: disable=broad-except,too-few-public-methods,line-too-long
class AnnoReaderBase(object):
    """Reads a Gene Association File. Returns a Python object."""

    tic = timeit.default_timer()

    # Expected values for a Qualifier
    exp_qualifiers = set([
        # Seen in both GAF and gene2go
        'not', 'contributes_to', 'colocalizes_with',
        # Seen in gene2go
        # TBD: resolve gaf(not|contributes_to) v gene2go(contributes_to)
        #'not contributes_to', 'not colocalizes_with',
        #
        # Although thee go not appear in:
        #     http://geneontology.org/page/go-annotation-conventions#qual
        # they do appear in more than one July 2018 GAFs:
        #     'enables', 'involved_in', 'part_of',
    ])

    def __init__(self, name, filename=None, **kws):
        # kws: allow_missing_symbol
        self.name = name
        self.filename = filename
        self.evobj = EvidenceCodes()
        # Read anotation file, store namedtuples:
        #     Gene2GoReader(filename=None, taxids=None):
        #     GafReader(filename=None, hdr_only=False, prt=sys.stdout, allow_missing_symbol=False):
        #     GpadReader(filename=None, hdr_only=False):
        self.hdr = None
        self.datobj = None
        # pylint: disable=no-member
        self.associations = self._init_associations(filename, **kws)
        # assert self.associations, 'NO ANNOTATIONS FOUND: {ANNO}'.format(ANNO=filename)

    def prt_summary_anno2ev(self, prt=sys.stdout):
        """Print annotation/evidence code summary."""
        self.evobj.prt_summary_anno2ev(self.associations, prt)

    def get_name(self):
        """Return type of annotation"""
        return self.name

    # pylint: disable=no-self-use
    def get_taxid(self):
        """Return taxid, if one was provided, otherwise return -1"""
        return -1

    def get_population(self):
        """Get population IDs (all DB_IDs)"""
        return self._get_population(self.associations)

    def get_ids_g_goids(self, goids):
        """Get database IDs (DB_IDs), given a set of GO IDs."""
        return set(nt.DB_ID for nt in self.associations if nt.GO_ID in goids)

    @staticmethod
    def _get_population(associations):
        """Get all IDs in the associations"""
        return set(nt.DB_ID for nt in associations)

    def get_id2gos(self, **kws):
        """Return all associations in a dict, id2gos"""
        return self._get_id2gos(self.associations, **kws)

    def _get_id2gos(self, associations, **kws):
        """Return given associations in a dict, id2gos"""
        options = AnnoOptions(self.evobj, **kws)
        # Default reduction is to remove. For all options, see goatools/anno/opts.py:
        #   * Evidence_Code == ND -> No biological data No biological Data available
        #   * Qualifiers contain NOT
        assc = self.reduce_annotations(associations, options)
        return self._get_dbid2goids(assc) if options.b_geneid2gos else self._get_goid2dbids(assc)

    # Qualifier (column 4)
    # Flags that modify the interpretation of an annotation one (or more) of NOT, contributes_to, colocalizes_with
    # This field is not mandatory;
    #     * cardinality 0, 1, >1;
    #     * for cardinality >1 use a pipe to separate entries (e.g. NOT|contributes_to)
    def prt_qualifiers(self, prt=sys.stdout):
        """Print Qualifiers: 1,462 colocalizes_with; 1,454 contributes_to; 1,157 not"""
        # 13 not colocalizes_with   (TBD: CHK - Seen in gene2go, but not gafs)
        #  4 not contributes_to     (TBD: CHK - Seen in gene2go, but not gafs)
        self._prt_qualifiers(self.associations, prt)

    @staticmethod
    def _prt_qualifiers(associations, prt=sys.stdout):
        """Print Qualifiers found in the annotations.
           QUALIFIERS:
                1,462 colocalizes_with
                1,454 contributes_to
                1,157 not
                   13 not colocalizes_with   (TBD: CHK - Seen in gene2go, but not gafs)
                    4 not contributes_to     (TBD: CHK - Seen in gene2go, but not gafs)
        """
        prt.write('QUALIFIERS:\n')
        for fld, cnt in cx.Counter(q for nt in associations for q in nt.Qualifier).most_common():
            prt.write('    {N:6,} {FLD}\n'.format(N=cnt, FLD=fld))

    def reduce_annotations(self, annotations, options):
        """Reduce annotations to ones used to identify enrichment (normally exclude ND and NOT)."""
        getfnc_qual_ev = options.getfnc_qual_ev()
        return [nt for nt in annotations if getfnc_qual_ev(nt.Qualifier, nt.Evidence_Code)]

    @staticmethod
    def _get_dbid2goids(associations):
        """Return gene2go data for user-specified taxids."""
        id2gos = cx.defaultdict(set)
        for ntd in associations:
            id2gos[ntd.DB_ID].add(ntd.GO_ID)
        return dict(id2gos)

    @staticmethod
    def _get_goid2dbids(associations):
        """Return gene2go data for user-specified taxids."""
        go2ids = cx.defaultdict(set)
        for ntd in associations:
            go2ids[ntd.GO_ID].add(ntd.DB_ID)
        return dict(go2ids)

    @staticmethod
    def get_date_yyyymmdd(yyyymmdd):
        """Return datetime.date given string."""
        return date(int(yyyymmdd[:4]), int(yyyymmdd[4:6], base=10), int(yyyymmdd[6:], base=10))

    def hms(self, msg, tic=None, prt=sys.stdout):
        """Print elapsed time and message."""
        if tic is None:
            tic = self.tic
        now = timeit.default_timer()
        hms = str(datetime.timedelta(seconds=(now-tic)))
        prt.write('{HMS}: {MSG}\n'.format(HMS=hms, MSG=msg))
        return now

    def chk_associations(self, fout_err=None):
        """Check that associations are in expected format."""
        pass

    def nts_ev_nd(self):
        """Get annotations where Evidence_code == 'ND' (No biological data)"""
        return [nt for nt in self.associations if nt.Evidence_Code == 'ND']

    def nts_qual_not(self):
        """Get annotations having Qualifiers containing NOT"""
        return [nt for nt in self.associations if self._has_not_qual(nt)]

    def chk_qualifiers(self):
        """Check format of qualifier"""
        if self.name == 'id2gos':
            return
        for ntd in self.associations:
            # print(ntd)
            qual = ntd.Qualifier
            assert isinstance(qual, set), '{NAME}: QUALIFIER MUST BE A LIST: {NT}'.format(
                NAME=self.name, NT=ntd)
            assert qual != set(['']), ntd
            assert qual != set(['-']), ntd
            assert 'always' not in qual, 'SPEC SAID IT WOULD BE THERE'

    @staticmethod
    def _has_not_qual(ntd):
        """Return True if the qualifiers contain a 'NOT'"""
        for qual in ntd.Qualifier:
            if 'not' in qual:
                return True
            if 'NOT' in qual:
                return True
        return False



# Copyright (C) 2016-2019, DV Klopfenstein, H Tang. All rights reserved."
