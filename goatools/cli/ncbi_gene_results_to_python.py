"""Read a NCBI Gene gene_result.txt file and write a Python module"""

from __future__ import print_function

__copyright__ = "Copyright (C) 2016-present, DV Klopfenstein, H Tang. All rights reserved."
__author__ = "DV Klopfenstein"

import os
from sys import stdout
import re
import datetime
import collections as cx
from argparse import ArgumentParser
from goatools.parsers.ncbi_gene_file_reader import NCBIgeneFileReader


def ncbi_tsv_to_py(fin_tsv, fout_py=None, prt=stdout):
    """Read a NCBI Gene file. Write data into one Python module per gene file"""
    obj = NCBIgeneToPythonCli()
    obj.ncbi_tsv_to_py(fin_tsv, fout_py, prt)

class NCBIgeneToPythonCli:
    """Read a NCBI Gene gene_result.txt file and write a Python module."""

    argparser = ArgumentParser(description='Convert a NCBI gene tsv file into a Python module')
    argparser.add_argument(
        'NCBI_gene_tsv', type=str, nargs='+',
        help='gene_result.tsv downloaded from NCBI Gene')
    argparser.add_argument(
        '-o', '--outfile',
        help='Write current citation report to an ASCII text file.')

    def cli(self, prt=stdout):
        """Command-line interface to print specified GO Terms from the DAG source ."""
        args = self.argparser.parse_args()
        # Aggregate all NCBI Gene data into a single output file
        if len(args.NCBI_gene_tsv) > 1 and args.outfile is not None:
            self.tsv_to_py_all(args.NCBI_gene_tsv, args.outfile, prt)
            return
        self.tsv_to_py_each(args.NCBI_gene_tsv, args.outfile, prt)

    def tsv_to_py(self, fin_tsv, fout_py=None, prt=stdout):
        """Read each NCBI Gene files. Write data into one Python module per gene file"""
        self.tsv_to_py_each([fin_tsv], fout_py, prt)

    def tsv_to_py_each(self, fin_tsvs, fout_py=None, prt=stdout):
        """Read each NCBI Gene files. Write data into one Python module per gene file"""
        in_outs = self._get_io_filenames(fin_tsvs, fout_py)
        for fin_tsv, fo_py in in_outs:
            self.ncbi_tsv_to_py(fin_tsv, fo_py, prt)

    def ncbi_tsv_to_py(self, fin_tsv, fout_py=None, prt=stdout):
        """Read a NCBI Gene file. Write data into one Python module per gene file"""
        nts = NCBIgeneFileReader(fin_tsv).get_nts()
        geneid2nt = self._get_geneid2nt(nts)
        self._wrpy_ncbi_gene_nts(fout_py, geneid2nt, prt)

    def _get_io_filenames(self, fin_tsvs, fout_py):
        """Get one output file for each input file"""
        nts = []
        ntobj = cx.namedtuple('NtIO', 'fin_tsv fout_py')
        ctr = cx.Counter()
        for fin_tsv in fin_tsvs:
            if os.path.exists(fin_tsv):
                if fout_py is not None:
                    nts.append(ntobj(fin_tsv=fin_tsv, fout_py=fout_py))
                else:
                    basename_tsv = os.path.basename(fin_tsv)
                    base, _ = os.path.splitext(basename_tsv)
                    fout_py_cur = self._get_foutpy(base, ctr[base])
                    ctr[base] += 1
                    nts.append(ntobj(fin_tsv=fin_tsv, fout_py=fout_py_cur))
            else:
                print('**ERROR-NOT FOUND: {FIN}'.format(FIN=fin_tsv))
        return nts

    @staticmethod
    def _get_foutpy(basename, cnt):
        """Get Python module name"""
        if cnt == 0:
            return '{F}.py'.format(F=basename)
        return '{F}{N}.py'.format(F=basename, N=cnt)

    def tsv_to_py_all(self, fin_tsvs, fout_py=None, prt=stdout):
        """Read all NCBI Gene files. Write all data into one Python module"""
        nts = self._read_tsvs_all(fin_tsvs)
        geneid2nt = self._get_geneid2nt(nts)
        self._wrpy_ncbi_gene_nts(fout_py, geneid2nt, prt)

    @staticmethod
    def _read_tsvs_all(fin_tsvs):
        """Read NCBI Gene tsv files. Return namedtuples"""
        nts_all = []
        for fin_tsv in fin_tsvs:
            if os.path.exists(fin_tsv):
                nts_cur = NCBIgeneFileReader(fin_tsv).get_nts()
                nts_all.extend(nts_cur)
            else:
                print('**ERROR-NOT FOUND: {FIN}'.format(FIN=fin_tsv))
        return nts_all

    @staticmethod
    def _get_geneid2nt(nts):
        """Get geneid2nt given a list of namedtuples."""
        geneid2nt = {}
        for ntd in nts:
            geneid = ntd.GeneID
            if geneid not in geneid2nt:
                geneid2nt[geneid] = ntd
            ## TBD: Some genes have more than one location
            ## elif self._nts_equal(geneid2nt[geneid], ntd):
            ##      print("DUPLICATE GeneID FOUND {N:9} {SYM}".format(N=geneid, SYM=ntd.Symbol))
        return geneid2nt

    @staticmethod
    def _nts_equal(nt0, nt1):
        """Return true if the namedtuples are equal"""
        return nt0 == nt1
        ## if nt0.Symbol[:3] != 'LOC':
        ##     if nt0 != nt1:
        ##         print('00000000000000000000', nt0)
        ##         print('11111111111111111111', nt1)
        ##     return nt0 == nt1
        ## else:
        ##     return True


    @staticmethod
    def _wrpy_ncbi_gene_nts(fout_py, geneid2nt, log):
        """Write namedtuples to a dict in a Python module."""
        num_genes = len(geneid2nt)
        with open(fout_py, 'w') as ofstrm:
            docstr = "Data downloaded from NCBI Gene converted into Python namedtuples."
            ofstrm.write('"""{PYDOC}"""\n\n'.format(PYDOC=docstr))
            ofstrm.write("from collections import namedtuple\n\n")
            ofstrm.write('WRITTEN = "{DATE}"'.format(
                DATE=re.sub('-', '_', str(datetime.date.today()))))
            ofstrm.write(' # {N} items\n\n'.format(N=num_genes))
            ntd = next(iter(geneid2nt.values())) # Access one dictionary value in Python 2
            ofstrm.write("#pylint: disable=line-too-long,too-many-lines,invalid-name\n")
            ofstrm.write("{NtName} = namedtuple('{NtName}', '{FLDS}')\n\n".format(
                NtName=type(ntd).__name__, FLDS=' '.join(ntd._fields)))
            ofstrm.write("GENEID2NT = {{ # {N:,} items\n".format(N=num_genes))
            for geneid, ntd in sorted(geneid2nt.items(), key=lambda t: t[0]):
                ofstrm.write("    {GeneID} : {NT},\n".format(GeneID=geneid, NT=ntd))
            ofstrm.write("}\n")
            log.write("  {N:10,} geneids WROTE: {PY}\n".format(N=num_genes, PY=fout_py))


# Copyright (C) 2016-present, DV Klopfenstein, H Tang. All rights reserved.
