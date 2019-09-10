"""Read a GO Annotation File (GAF) and store the data in a Python object.

    Annotations available from the Gene Ontology Consortium:
        http://current.geneontology.org/annotations/
"""

import sys
from goatools.anno.annoreader_base import AnnoReaderBase
from goatools.anno.init.reader_gaf import GafData
from goatools.anno.init.reader_gaf import InitAssc

__copyright__ = "Copyright (C) 2016-2019, DV Klopfenstein, H Tang. All rights reserved."
__author__ = "DV Klopfenstein"


# pylint: disable=broad-except,too-few-public-methods,line-too-long
class GafReader(AnnoReaderBase):
    """Reads a Gene Annotation File (GAF). Returns a Python object."""

    exp_kws = {'hdr_only', 'prt', 'namespaces', 'allow_missing_symbol', 'godag'}

    def __init__(self, filename=None, **kws):
        super(GafReader, self).__init__(
            'gaf', filename,
            godag=kws.get('godag'),
            hdr_only=kws.get('hdr_only', False),
            prt=kws.get('prt', sys.stdout),
            namespaces=kws.get('namespaces'),
            allow_missing_symbol=kws.get('allow_missing_symbol', False))

    def read_gaf(self, namespace='BP', **kws):
        """Read Gene Association File (GAF). Return associations."""
        return self.get_id2gos(namespace, **kws)

    @staticmethod
    def wr_txt(fout_gaf, nts):
        """Write namedtuples into a gaf format"""
        pat = (
            '{DB}\t{DB_ID}\t{DB_Symbol}\t{Qualifier}\t{GO_ID}\t{DB_Reference}\t'
            '{Evidence_Code}\t{With_From}\t{NS}\t{DB_Name}\t{DB_Synonym}\t{DB_Type}\t'
            '{Taxon}\t{Date}\t{Assigned_By}\t{Extension}\t{Gene_Product_Form_ID}\n')
        sets = {'Qualifier', 'DB_Reference', 'With_From', 'DB_Name', 'DB_Synonym', 'Gene_Product_Form_ID'}
        ns2a = {ns:p for p, ns in GafData.aspect2ns.items()}
        with open(fout_gaf, 'w') as prt:
            prt.write('!gaf-version: 2.1\n')
            for ntd in nts:
                dct = ntd._asdict()
                for fld in sets:
                    dct[fld] = '|'.join(sorted(dct[fld]))
                dct['Taxon'] = '|'.join(['taxon:{T}'.format(T=t) for t in dct['Taxon']])
                dct['NS'] = ns2a[dct['NS']]
                dct['Date'] = dct['Date'].strftime('%Y%m%d')
                prt.write(pat.format(**dct))
                #prt.write('{NT}\n'.format(NT=ntd))
        print('  {N} annotations WROTE: {GAF}'.format(N=len(nts), GAF=fout_gaf))

    def chk_associations(self, fout_err="gaf.err"):
        """Check that fields are legal in GAF"""
        obj = GafData("2.1")
        return obj.chk(self.associations, fout_err)

    def has_ns(self):
        """Return True if namespace field, NS exists on annotation namedtuples"""
        return True

    # def _init_associations(self, fin_gaf, hdr_only, prt, namespaces, allow_missing_symbol):
    def _init_associations(self, fin_gaf, **kws):
        """Read annotation file and store a list of namedtuples."""
        ini = InitAssc(fin_gaf)
        nts = ini.init_associations(kws['hdr_only'], kws['prt'], kws['namespaces'], kws['allow_missing_symbol'])
        self.hdr = ini.hdr
        return nts


# Copyright (C) 2016-2019, DV Klopfenstein, H Tang. All rights reserved."
