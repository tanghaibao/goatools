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

__copyright__ = "Copyright (C) 2016-2019, DV Klopfenstein, H Tang. All rights reserved."
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

    def get_ns2ntsanno(self, taxid=None):
        """Return all associations in three (one for BP MF CC) dicts, id2gos"""
        # kws1: taxid
        ntsanno = self.get_associations(taxid)
        # kws2: ev_include ev_exclude ...
        return self._get_ns2ntsanno(ntsanno)

    def get_associations(self, taxid=None):
        """Return annotations"""
        # kws: taxid
        if len(self.taxid2asscs) == 1:
            taxid_cur = next(iter(self.taxid2asscs.keys()))
            return self.taxid2asscs[taxid_cur]
        assert taxid in self.taxid2asscs, '**FATAL: TAXID({T}) DATA MISSING'.format(T=taxid)
        return self.taxid2asscs[taxid]

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
    @staticmethod
    def _init_associations(fin_anno, taxid=None, taxids=None, namespaces=None, **kws):
        """Read annotation file and store a list of namedtuples."""
        return InitAssc(taxid, taxids).init_associations(fin_anno, taxids, namespaces)

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


# Copyright (C) 2016-2019, DV Klopfenstein, H Tang. All rights reserved."
