"""Read a NCBI Gene gene_result.txt file and write a Python module.

Usage:
  ncbi_gene_results_to_python.py [options]

Options:
  -h --help                                 show this help message and exit

  -i <gene_result.txt>   Read NCBI Gene file [default: gene_result.txt]
  -o <gene_result.py>    Write Python file [default: gene_result.py]

"""

from __future__ import print_function

__copyright__ = "Copyright (C) 2016-2018, DV Klopfenstein, H Tang. All rights reserved."
__author__ = "DV Klopfenstein"

import os
import sys
import re
import datetime
from goatools.cli.docopt_parse import DocOptParse
from goatools.parsers.ncbi_gene_file_reader import NCBIgeneFileReader


# pylint: disable=too-few-public-methods
class NCBIgeneToPythonCli(object):
    """Read a NCBI Gene gene_result.txt file and write a Python module."""

    kws_dict = set(['i', 'o'])

    def __init__(self):
        self.objdoc = DocOptParse(__doc__, self.kws_dict, set())

    def cli(self, prt=sys.stdout):
        """Command-line interface to print specified GO Terms from the DAG source ."""
        kws = self.objdoc.get_docargs(prt=None)
        if os.path.exists(kws['i']):
            obj = NCBIgeneFileReader(kws['i'])
            nts = obj.get_nts()
            if nts:
                geneid2nt = self._get_geneid2nt(nts)
                self._wrpy_ncbi_gene_nts(kws['o'], geneid2nt, prt)
        else:
            raise RuntimeError("\n{DOC}\n**ERROR: NO FILE FOUND: {NCBI}".format(
                NCBI=kws['i'], DOC=__doc__))

    @staticmethod
    def _get_geneid2nt(nts):
        """Get geneid2nt given a list of namedtuples."""
        geneid2nt = {}
        for ntd in nts:
            geneid = ntd.GeneID
            if geneid not in geneid2nt:
                geneid2nt[geneid] = ntd
            else:
                print("DUPLICATE GeneID FOUND {N:9} {SYM}".format(N=geneid, SYM=ntd.Symbol))
        return geneid2nt

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
            log.write("  {N:9} geneids WROTE: {PY}\n".format(N=num_genes, PY=fout_py))


# Copyright (C) 2016-2018, DV Klopfenstein, H Tang. All rights reserved.
