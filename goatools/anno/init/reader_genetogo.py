"""Read an NCBI gene2go GO Association File and store the data in a Python object.

    Annotations available from NCBI:
        ftp://ftp.ncbi.nih.gov/gene/DATA/gene2go.gz

"""

import sys
import collections as cx
import timeit
import datetime

__copyright__ = "Copyright (C) 2016-2019, DV Klopfenstein, H Tang. All rights reserved."
__author__ = "DV Klopfenstein"


# pylint: disable=too-few-public-methods
class InitAssc(object):
    """Read annotation file and store a list of namedtuples."""

    taxid_dflt = 9606  # human is the default taxid

    # Shared names(3): GO_ID Qualifier GO_term
    # Equivalent names(3): DB_ID->GeneID  Evidence_Code->Evidence  DB_Reference->PubMed
    # Only gene2go(2): Category/NS and tax_id
    # pylint: disable=line-too-long
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
    def init_associations(self, fin_anno, taxids=None, namespaces=None):
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
                get_all_taxids = taxids is True
                get_all_nss = namespaces is None or namespaces == {'BP', 'MF', 'CC'}
                taxids = self.taxids
                for lnum, line in enumerate(ifstrm, 1):
                    # Read data
                    if line[0] != '#':
                        vals = line.split('\t')
                        taxid = int(vals[0])
                        nspc = category2ns[vals[7].rstrip()]
                        if (get_all_taxids or taxid in taxids) and (get_all_nss or nspc in namespaces):
                            # assert len(vals) == 8
                            ntd = ntobj(
                                tax_id=taxid,
                                DB_ID=int(vals[1]),
                                GO_ID=vals[2],
                                Evidence_Code=vals[3],
                                Qualifier=self._get_qualifiers(vals[4]),
                                GO_term=vals[5],
                                DB_Reference=self._get_pmids(vals[6]),
                                NS=nspc)
                            #self._chk_qualifiers(qualifiers, lnum, ntd)
                            nts.append(ntd)
                    # Read header
                    elif line[0] == '#':
                        assert line[1:-1].split('\t') == self.hdrs
        # pylint: disable=broad-except
        except Exception as inst:
            import traceback
            traceback.print_exc()
            sys.stderr.write("\n  **FATAL: {MSG}\n\n".format(MSG=str(inst)))
            sys.stderr.write("**FATAL: {FIN}[{LNUM}]:\n{L}".format(FIN=fin_anno, L=line, LNUM=lnum))
            self._prt_line_detail(sys.stdout, line, lnum)
            sys.exit(1)
        print('HMS:{HMS} {N:7,} annotations READ: {ANNO} {NSs}'.format(
            N=len(nts), ANNO=fin_anno,
            NSs=','.join(namespaces) if namespaces else '',
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
