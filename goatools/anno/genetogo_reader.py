"""Read an NCBI gene2go GO Association File and store the data in a Python object.

    Annotations available from NCBI:
        ftp://ftp.ncbi.nih.gov/gene/DATA/gene2go.gz

"""

import sys
import collections as cx
import timeit
import datetime
from goatools.anno.annoreader_base import AnnoReaderBase

__copyright__ = "Copyright (C) 2016-2019, DV Klopfenstein, H Tang. All rights reserved."
__author__ = "DV Klopfenstein"


# pylint: disable=broad-except,too-few-public-methods,line-too-long
class Gene2GoReader(AnnoReaderBase):
    """Reads a Gene Annotation File (GAF). Returns a Python object."""

    taxid_dflt = 9606  # human is the default taxid
    # Shared names(3): GO_ID Qualifier GO_term
    # Equivalent names(3): DB_ID->GeneID  Evidence_Code->Evidence  DB_Reference->PubMed
    # Only gene2go(2): Category/NS and tax_id
    hdrs = ['tax_id', 'GeneID', 'GO_ID', 'Evidence', 'Qualifier', 'GO_term', 'PubMed', 'Category']
    flds = ['tax_id', 'DB_ID', 'GO_ID', 'Evidence_Code', 'Qualifier', 'GO_term', 'DB_Reference', 'NS']

    ## exp_kwdct = set(['allow_missing_symbol'])

    def __init__(self, filename=None, taxids=None, prt=sys.stdout):  # , **kws):
        super(Gene2GoReader, self).__init__(filename, prt)
        # kws: allow_missing_symbol
        ## self.kws = {k:v for k, v in kws.items() if k in self.exp_kwdct}
        #### self.filename = filename
        #### self.evobj = EvidenceCodes()
        # Initialize associations and header information
        #### self.hdr = None
        #### self.datobj = None
        self.associations = self._read_nts(filename, taxids)
        self.taxid2asscs = self._init_taxid2asscs(self.associations)
        #### self.associations = self.evobj.sort_nts(self._read_nts(filename), 'Evidence_Code')
        #### self.associations = self._init_assn(filename, hdr_only, prt) if filename is not None else []

    def prt_qualifiers(self, prt=sys.stdout):
        """Print Qualifiers found in the annotations.
           QUALIFIERS:
                1,462 colocalizes_with
                1,454 contributes_to
                1,157 not
                   13 not colocalizes_with   (TBD: CHK - Seen in gene2go, but not gafs)
                    4 not contributes_to     (TBD: CHK - Seen in gene2go, but not gafs)
        """
        prt.write('QUALIFIERS:\n')
        for fld, cnt in cx.Counter(q for nt in self.associations for q in nt.Qualifier).most_common():
            prt.write('    {N:6,} {FLD}\n'.format(N=cnt, FLD=fld))

    def get_annotations_dct(self, taxid, options):
        """Return gene2go data for user-specified taxids."""
        # Simple associations
        id2gos = cx.defaultdict(set)
        b_geneid2gos = options.b_geneid2gos
        keep = options.keep
        assert taxid in self.taxid2asscs, '**FATAL: TAXID({T}) DATA MISSING'.format(T=taxid)
        for ntd in self.taxid2asscs[taxid]:
            # NOT: Used when gene is expected to have function F, but does NOT.
            # ND : GO function not seen after exhaustive annotation attempts to the gene.
            # if 'not' not in set(ntd.Qualifier) and ntd.Evidence_Code != 'ND':
            if keep(ntd.Qualifier, ntd.Evidence_Code):
                if b_geneid2gos:
                    id2gos[int(ntd.DB_ID)].add(ntd.GO_ID)
                else:
                    id2gos[ntd.GO_ID].add(int(ntd.DB_ID))
        return dict(id2gos)

    def get_annotations_taxid2dct(self, options, taxids=None):
        """Read Gene Association File (GAF). Return data."""
        taxids_stored = set(self.taxid2asscs.keys())
        if taxids is None:
            taxids = taxids_stored
        else:
            taxids = set(taxids).intersection(taxids_stored)
        keep = options.keep
        taxid2asscs = cx.defaultdict(lambda: cx.defaultdict(lambda: cx.defaultdict(set)))
        for taxid in taxids:
            for ntd in self.taxid2asscs[taxid]:
                if keep(ntd.Qualifier, ntd.Evidence_Code):
                    geneid = ntd.DB_ID
                    go_id = ntd.GO_ID
                    if taxid:
                        taxid2asscs[taxid]['GeneID2GOs'][geneid].add(go_id)
                        taxid2asscs[taxid]['GO2GeneIDs'][go_id].add(geneid)
        return taxid2asscs

    @staticmethod
    def fill_taxid2asscs(taxid2asscs_usr, taxid2asscs_ret):
        """Fill user taxid2asscs for backward compatibility."""
        for taxid, ab_ret in taxid2asscs_ret.items():
            taxid2asscs_usr[taxid]['GeneID2GOs'] = ab_ret['GeneID2GOs']
            taxid2asscs_usr[taxid]['GO2GeneIDs'] = ab_ret['GO2GeneIDs']

    @staticmethod
    def get_id2gos_all(taxid2asscs_a2b):
        """Get associations for all stored species taxid2asscs[taxid][GeneID2GOs|GO2GeneIDs]."""
        id2gos_all = {}
        for a2b in taxid2asscs_a2b.values():
            for geneid, gos in a2b['GeneID2GOs'].items():
                id2gos_all[geneid] = gos
        return id2gos_all

    ## def _init_assn(self, fin_anno, hdr_only, prt):
    ##     """Read GAF file. Store annotation data in a list of namedtuples."""
    ##     nts = self._read_nts(fin_anno, hdr_only)
    ##     # Annotation file has been read
    ##     if prt:
    ##         prt.write("  READ    {N:9,} associations: {FIN}\n".format(N=len(nts), FIN=fin_anno))
    ##     # If there are illegal GAF lines ...
    ##     if self.datobj:
    ##         if self.datobj.ignored or self.datobj.illegal_lines:
    ##             self.datobj.prt_error_summary(fin_anno)
    ##     return self.evobj.sort_nts(nts, 'Evidence_Code')

    def _chk_qualifiers(self, qualifiers, lnum, ntd):
        """Check that qualifiers are expected values."""
        # http://geneontology.org/page/go-annotation-conventions#qual
        for qualifier in qualifiers:
            if qualifier not in self.exp_qualifiers:
                errname = '** WARNING: UNEXPECTED QUALIFIER({QUAL})'.format(QUAL=qualifier)
                print('LNUM({LNUM}): {ERR}\n{NT}'.format(LNUM=lnum, ERR=errname, NT=ntd))
                # self.illegal_lines[errname].append((lnum, "\t".join(flds)))

    def _init_taxid2asscs(self, associations):
        """Create dict with taxid keys and annotation namedtuple list."""
        taxid2asscs = cx.defaultdict(list)
        # Loop through list of a namedtuples saved from NCBI's gene2go annotations
        for ntanno in associations:
            taxid2asscs[ntanno.tax_id].append(ntanno)
        print('{N} taxids stored'.format(N=len(taxid2asscs)))
        return {t:self.evobj.sort_nts(nts, 'Evidence_Code') for t, nts in taxid2asscs.items()}

    def _read_nts(self, fin_anno, taxids=None):
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
                taxids = self._get_taxids(taxids, sys.stdout)
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
        print('HMS:{HMS} {N:,} lines READ: {ANNO}'.format(
            N=len(nts), ANNO=fin_anno,
            HMS=str(datetime.timedelta(seconds=(timeit.default_timer()-tic)))))
        return nts

    def _get_taxids(self, taxids, prt=None):
        """Get a function which determines whether to save annotations for a given taxid."""
        # Default taxid is Human
        if taxids is None:
            if prt:
                prt.write('**NOTE: DEFAULT TAXID STORED FROM gene2go IS 9606 (human)\n')
            return set([self.taxid_dflt])
        if isinstance(taxids, int):
            return set([taxids])
        assert hasattr(taxids, '__iter__'), 'BAD TAXIDS({T})'.format(T=taxids)
        return set(taxids)

    @staticmethod
    def _get_qualifiers(qualifier):
        """Return a list of qualifiers if they exist."""
        if qualifier == '-':
            return {}
        # TBD: Separator '|' or ' '
        return {q.lower() for q in qualifier.split('|')}

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


# Copyright (C) 2016-2019, DV Klopfenstein, H Tang. All rights reserved."
