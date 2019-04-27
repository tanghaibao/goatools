"""Read an NCBI gene2go GO Association File and store the data in a Python object.

    Annotations available from NCBI:
        ftp://ftp.ncbi.nih.gov/gene/DATA/gene2go.gz

"""

import sys
import collections as cx
import timeit
import datetime
from goatools.anno.annoreader_base import AnnoReaderBase
from goatools.anno.opts import AnnoOptions

__copyright__ = "Copyright (C) 2016-2019, DV Klopfenstein, H Tang. All rights reserved."
__author__ = "DV Klopfenstein"


# pylint: disable=broad-except,too-few-public-methods,line-too-long
class Gene2GoReader(AnnoReaderBase):
    """Reads a Gene Annotation File (GAF). Returns a Python object."""

    def __init__(self, filename=None, **kws):
        # kws: taxids or taxid
        super(Gene2GoReader, self).__init__('gene2go', filename, **kws)
        # Initialize associations and header information
        self.taxid2asscs = self._init_taxid2asscs()

    def get_id2gos(self, **kws):
    #### def get_annotations_dct(self, taxid, options):
        """Return geneid2gos, or optionally go2geneids."""
        if len(self.taxid2asscs) == 1:
            taxid = next(iter(self.taxid2asscs.keys()))
            return self._get_id2gos(self.taxid2asscs[taxid], **kws)
        assert 'taxid' in kws, "**FATAL: 'taxid' NOT FOUND IN Gene2GoReader::get_id2gos({KW})".format(KW=kws)
        taxid = kws['taxid']
        assert taxid in self.taxid2asscs, '**FATAL: TAXID({T}) DATA MISSING'.format(T=taxid)
        return self._get_id2gos(self.taxid2asscs[taxid], **kws)

    def get_name(self):
        """Get name using taxid"""
        if len(self.taxid2asscs) == 1:
            return '{BASE}_{TAXID}'.format(
                BASE=self.name, TAXID=next(iter(self.taxid2asscs.keys())))
        return '{BASE}_various'.format(BASE=self.name)

    def get_taxid(self):
        """Return taxid, if one was provided. Other wise return True representing all taxids"""
        return next(iter(self.taxid2asscs.keys())) if len(self.taxid2asscs) == 1 else True

    # -- taxids2asscs -------------------------------------------------------------------------
    def get_taxid2asscs(self, taxids=None, **kws):
        """Read Gene Association File (GAF). Return data."""
        # WAS: get_annotations_taxid2dct
        taxid2asscs = cx.defaultdict(lambda: cx.defaultdict(lambda: cx.defaultdict(set)))
        options = AnnoOptions(self.evobj, **kws)
        for taxid in self._get_taxids(taxids):
            nts = self.taxid2asscs[taxid]
            assc = self.reduce_annotations(nts, options)
            taxid2asscs[taxid]['ID2GOs'] = self._get_dbid2goids(assc)
            taxid2asscs[taxid]['GO2IDs'] = self._get_goid2dbids(assc)
        return taxid2asscs

    @staticmethod
    def fill_taxid2asscs(taxid2asscs_usr, taxid2asscs_ret):
        """Fill user taxid2asscs for backward compatibility."""
        for taxid, ab_ret in taxid2asscs_ret.items():
            taxid2asscs_usr[taxid]['ID2GOs'] = ab_ret['ID2GOs']
            taxid2asscs_usr[taxid]['GO2IDs'] = ab_ret['GO2IDs']

    @staticmethod
    def get_id2gos_all(taxid2asscs_a2b):
        """Get associations for all stored species taxid2asscs[taxid][ID2GOs|GO2IDs]."""
        id2gos_all = {}
        for a2b in taxid2asscs_a2b.values():
            for geneid, gos in a2b['ID2GOs'].items():
                id2gos_all[geneid] = gos
        return id2gos_all

    def _get_taxids(self, taxids=None):
        """Return user-specified taxids or taxids in self.taxid2asscs"""
        taxid_keys = set(self.taxid2asscs.keys())
        return taxid_keys if taxids is None else set(taxids).intersection(taxid_keys)

    # -- initialization -----------------------------------------------------------------------
    @staticmethod
    def _init_associations(fin_anno, taxid=None, taxids=None):
        """Read annotation file and store a list of namedtuples."""
        return _InitAssc(taxid, taxids).init_associations(fin_anno, taxids)

    def _init_taxid2asscs(self):
        """Create dict with taxid keys and annotation namedtuple list."""
        taxid2asscs = cx.defaultdict(list)
        for ntanno in self.associations:
            taxid2asscs[ntanno.tax_id].append(ntanno)
        assert len(taxid2asscs) != 0, "**FATAL: NO TAXIDS: {F}".format(F=self.filename)
        # """Print the number of taxids stored."""
        prt = sys.stdout
        num_taxids = len(taxid2asscs)
        prt.write('{N} taxids stored'.format(N=num_taxids))
        if num_taxids < 5:
            prt.write(': {Ts}'.format(Ts=' '.join(sorted(str(t) for t in taxid2asscs))))
        prt.write('\n')
        return dict(taxid2asscs)


class _InitAssc(object):
    """Read annotation file and store a list of namedtuples."""

    taxid_dflt = 9606  # human is the default taxid

    # Shared names(3): GO_ID Qualifier GO_term
    # Equivalent names(3): DB_ID->GeneID  Evidence_Code->Evidence  DB_Reference->PubMed
    # Only gene2go(2): Category/NS and tax_id
    hdrs = ['tax_id', 'GeneID', 'GO_ID', 'Evidence', 'Qualifier', 'GO_term', 'PubMed', 'Category']
    flds = ['tax_id', 'DB_ID', 'GO_ID', 'Evidence_Code', 'Qualifier', 'GO_term', 'DB_Reference', 'NS']

    def __init__(self, taxid=None, taxids=None):
        self.taxids = self._init_taxids(taxid, taxids)

    @staticmethod
    def _init_taxids(taxid, taxids):
        """Return taxid set"""
        ret = set()
        if taxids is not None:
            if taxids is True:
                return True
            if isinstance(taxids, int):
                ret.add(taxids)
            else:
                ret.update(taxids)
        if taxid is not None:
            ret.add(taxid)
        if not ret:
            ret.add(9606)
            # pylint: disable=superfluous-parens
            print('**NOTE: DEFAULT TAXID STORED FROM gene2go IS 9606 (human)\n')
        return ret

    # pylint: disable=too-many-locals
    def init_associations(self, fin_anno, taxids=None):
        """Read annotation file. Store annotation data in a list of namedtuples."""
        nts = []
        if fin_anno is None:
            return nts
        tic = timeit.default_timer()
        lnum = -1
        line = "\t"*len(self.flds)
        try:
            with open(fin_anno) as ifstrm:
                category2ns = {'Process':'BP', 'Function':'MF', 'Component':'CC'}
                ntobj = cx.namedtuple('ntanno', self.flds)
                # Get: 1) Specified taxids, default taxid(human), or all taxids
                get_all = taxids is True
                taxids = self.taxids
                for lnum, line in enumerate(ifstrm, 1):
                    # Read data
                    if line[0] != '#':
                        vals = line.split('\t')
                        taxid = int(vals[0])
                        if get_all or taxid in taxids:
                            # assert len(vals) == 8
                            ntd = ntobj(
                                tax_id=taxid,
                                DB_ID=int(vals[1]),
                                GO_ID=vals[2],
                                Evidence_Code=vals[3],
                                Qualifier=self._get_qualifiers(vals[4]),
                                GO_term=vals[5],
                                DB_Reference=self._get_pmids(vals[6]),
                                NS=category2ns[vals[7].rstrip()])
                            #self._chk_qualifiers(qualifiers, lnum, ntd)
                            nts.append(ntd)
                    # Read header
                    elif line[0] == '#':
                        assert line[1:-1].split('\t') == self.hdrs
        except Exception as inst:
            import traceback
            traceback.print_exc()
            sys.stderr.write("\n  **FATAL: {MSG}\n\n".format(MSG=str(inst)))
            sys.stderr.write("**FATAL: {FIN}[{LNUM}]:\n{L}".format(FIN=fin_anno, L=line, LNUM=lnum))
            self._prt_line_detail(sys.stdout, line, lnum)
            sys.exit(1)
        print('HMS:{HMS} {N:7,} annotations READ: {ANNO}'.format(
            N=len(nts), ANNO=fin_anno,
            HMS=str(datetime.timedelta(seconds=(timeit.default_timer()-tic)))))
        return nts

    @staticmethod
    def _get_qualifiers(qualifier):
        """Return a list of qualifiers if they exist."""
        return set(qualifier.split(' ')) if qualifier != '-' else set()

    @staticmethod
    def _get_pmids(pmidstr):
        """Return a list of PMIDs if they exist."""
        if pmidstr == '-':
            return []
        return ['PMID:{N}'.format(N=n) for n in pmidstr.split('|')]

    def _prt_line_detail(self, prt, line, lnum=""):
        """Print each field and its value."""
        data = zip(self.flds, line.split('\t'))
        txt = ["{:2}) {:13} {}".format(i, hdr, val) for i, (hdr, val) in enumerate(data)]
        prt.write("{LNUM}\n{TXT}\n".format(LNUM=lnum, TXT='\n'.join(txt)))

    ## def _chk_qualifiers(self, qualifiers, lnum, ntd):
    ##     """Check that qualifiers are expected values."""
    ##     # http://geneontology.org/page/go-annotation-conventions#qual
    ##     for qualifier in qualifiers:
    ##         if qualifier not in self.exp_qualifiers:
    ##             errname = '** WARNING: UNEXPECTED QUALIFIER({QUAL})'.format(QUAL=qualifier)
    ##             # pylint: disable=superfluous-parens
    ##             print('LNUM({LNUM}): {ERR}\n{NT}'.format(LNUM=lnum, ERR=errname, NT=ntd))
    ##             # self.illegal_lines[errname].append((lnum, "\t".join(flds)))


# Copyright (C) 2016-2019, DV Klopfenstein, H Tang. All rights reserved."
