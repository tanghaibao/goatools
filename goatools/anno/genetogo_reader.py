"""Read an NCBI gene2go GO Association File and store the data in a Python object.

    Annotations available from NCBI:
        ftp://ftp.ncbi.nih.gov/gene/DATA/gene2go.gz

"""

import sys
import collections as cx
from itertools import chain
from goatools.anno.init.reader_genetogo import InitAssc
from goatools.anno.annoreader_base import AnnoReaderBase
from goatools.anno.opts import AnnoOptions

__copyright__ = "Copyright (C) 2016-present, DV Klopfenstein, H Tang. All rights reserved."
__author__ = "DV Klopfenstein"


# pylint: disable=broad-except,too-few-public-methods,line-too-long
class Gene2GoReader(AnnoReaderBase):
    """Reads a Gene Annotation File (GAF). Returns a Python object."""

    exp_kws = {'taxids', 'taxid', 'namespaces', 'godag'}

    def __init__(self, filename=None, **kws):
        # kws: taxids or taxid
        super(Gene2GoReader, self).__init__('gene2go', filename, **kws)
        # Each taxid has a list of namedtuples - one for each line in the annotations
        self.taxid2asscs = self._init_taxid2asscs()

    def get_ns2assc(self, taxid=None, **kws):
        """Return given associations into 3 (BP, MF, CC) dicts, id2gos"""
        return {ns:self._get_id2gos(nts, **kws) for ns, nts in self.get_ns2ntsanno(taxid).items()}

    def get_ns2ntsanno(self, taxid=None):
        """Return all associations in three (one for BP MF CC) dicts, id2gos"""
        ntsanno = self.get_associations(taxid)
        # kws2: ev_include ev_exclude ...
        return self._get_ns2ntsanno(ntsanno)

    def get_associations(self, taxid=None):
        """Return annotations"""
        # If only one taxid is loaded, ignore the taxid arg
        if len(self.taxid2asscs) == 1:
            taxid_cur = next(iter(self.taxid2asscs.keys()))
            return self.taxid2asscs[taxid_cur]
        # If taxid is True, combine all associations
        if taxid is True:
            return list(chain.from_iterable(self.taxid2asscs.values()))
        # If taxid is an int, return the associations for user-specified taxid
        if isinstance(taxid, int):
            return self.taxid2asscs[taxid] if taxid in self.taxid2asscs else []
        # If multiple taxids were loaded, taxid, must be specified
        if taxid is None or not taxid:
            return self._warning_taxid(taxid)
        # Assume taxid is a list or set and combine those taxids
        taxids = set(taxid).intersection(self.taxid2asscs.keys())
        if taxids:
            # Return user-specified taxids combined
            return list(chain.from_iterable(self.taxid2asscs[t] for t in taxids))
        return {}

    @staticmethod
    def _warning_taxid(taxid):
        """Warn if an unexpected taxid"""
        pat = ('**WARNING: NO ASSOCIATIONS FOR taxid({TAXID}). '
               'Taxid MUST BE AN int, list of ints, OR bool')
        print(pat.format(TAXID=taxid))
        return {}

    def get_id2gos_nss(self, **kws):
        """Return all associations in a dict, id2gos, regardless of namespace"""
        taxids = self._get_taxids(kws.get('taxids'), kws.get('taxid'))
        assert taxids, "NO TAXIDS FOUND"
        assc = list(chain.from_iterable([self.taxid2asscs[t] for t in taxids]))
        return self._get_id2gos(assc, **kws)

    def get_name(self):
        """Get name using taxid"""
        if len(self.taxid2asscs) == 1:
            return '{BASE}_{TAXID}'.format(
                BASE=self.name, TAXID=next(iter(self.taxid2asscs.keys())))
        return '{BASE}_various'.format(BASE=self.name)

    def get_taxid(self):
        """Return taxid, if one was provided. Other wise return True representing all taxids"""
        return next(iter(self.taxid2asscs.keys())) if len(self.taxid2asscs) == 1 else True

    def has_ns(self):
        """Return True if namespace field, NS exists on annotation namedtuples"""
        return True

    def prt_counts(self, prt=sys.stdout):
        """Print the number of taxids stored."""
        num_taxids = len(self.taxid2asscs)
        num_annos = sum(len(a) for a in self.taxid2asscs.values())
        # 792,891 annotations for 3 taxids stored: 10090 7227 9606
        cnts = self._get_counts(list(chain.from_iterable(self.taxid2asscs.values())))
        prt.write('{A:8,} annotations, {P:,} proteins/genes, {G:,} GO IDs, {N} taxids stored'.format(
            A=num_annos, N=num_taxids, G=cnts['GOs'], P=cnts['geneids']))
        if num_taxids < 5:
            prt.write(': {Ts}'.format(Ts=' '.join(str(t) for t in sorted(self.taxid2asscs))))
        prt.write('\n')
        # 102,430 annotations for taxid  7227
        # 323,776 annotations for taxid  9606
        # 366,685 annotations for taxid 10090
        if num_taxids == 1:
            return
        for taxid, assc in self.taxid2asscs.items():
            cnts = self._get_counts(assc)
            prt.write('{A:8,} annotations, {P:,} proteins/genes, {G:,} GO IDs for taxid {T}\n'.format(
                A=len(assc), T=taxid, G=cnts['GOs'], P=cnts['geneids']))

    @staticmethod
    def _get_counts(nts):
        """Return the count of GO IDs and genes/proteins in a set of annotation namedtuples"""
        sets = cx.defaultdict(set)
        for ntd in nts:
            sets['geneids'].add(ntd.DB_ID)
            sets['GOs'].add(ntd.GO_ID)
        return {'GOs':len(sets['GOs']), 'geneids':len(sets['geneids'])}

    # -- taxids2asscs -------------------------------------------------------------------------
    def get_taxid2asscs(self, taxids=None, **kws):
        """Read Gene Association File (GAF). Return data."""
        # WAS: get_annotations_taxid2dct
        taxid2asscs = cx.defaultdict(lambda: cx.defaultdict(lambda: cx.defaultdict(set)))
        options = AnnoOptions(self.evobj, **kws)
        for taxid in self._get_taxids(taxids):
            nts = self.taxid2asscs[taxid]
            assc = self.reduce_annotations(nts, options)
            taxid2asscs[taxid]['ID2GOs'] = self.get_dbid2goids(assc)
            taxid2asscs[taxid]['GO2IDs'] = self.get_goid2dbids(assc)
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

    def _get_taxids(self, taxids=None, taxid=None):
        """Return user-specified taxids or taxids in self.taxid2asscs"""
        taxid_keys = set(self.taxid2asscs.keys())
        if taxids is None and taxid is not None:
            taxids = [taxid]
        return taxid_keys if taxids is None else set(taxids).intersection(taxid_keys)

    # -- initialization -----------------------------------------------------------------------
    # pylint: disable=unused-argument
    @staticmethod
    def _init_associations(fin_anno, taxid=None, taxids=None, namespaces=None, **kws):
        """Read annotation file and store a list of namedtuples."""
        return InitAssc(taxid, taxids).init_associations(fin_anno, taxids, namespaces)

    def _init_taxid2asscs(self):
        """Create dict with taxid keys and annotation namedtuple list."""
        taxid2asscs = cx.defaultdict(list)
        for ntanno in self.associations:
            taxid2asscs[ntanno.tax_id].append(ntanno)
        assert taxid2asscs, "**FATAL: NO TAXIDS: {F}".format(F=self.filename)
        return dict(taxid2asscs)


# Copyright (C) 2016-present, DV Klopfenstein, H Tang. All rights reserved."
