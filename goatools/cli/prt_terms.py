"""Command-line interface to print specified GO Terms from the DAG source

Usage:
  prt_terms.py [GO ...] [GO_FILE]
  prt_terms.py [GO ...] [GO_FILE] [options]

Options:
  -h --help                                 show this help message and exit

  -i <file.txt>, --ifile=<sections_in.txt>  Read or Write file name [default: sections_in.txt]
  -n <GO_name>, --name=<GO_name>            Name of a GO Term

  --obo=<file.obo>     Ontologies in obo file [default: go-basic.obo].
"""

from __future__ import print_function

__copyright__ = "Copyright (C) 2016-2018, DV Klopfenstein, H Tang. All rights reserved."
__author__ = "DV Klopfenstein"


import sys
from goatools.cli.docopt_parse import DocOptParse
from goatools.cli.gos_get import GetGOs
from goatools.test_data.wr_subobo import WrSubObo


# pylint: disable=too-few-public-methods
class PrtGOterms(object):
    """Command-line interface to print specified GO Terms from the DAG source."""

    kws_dict = set(['GO', 'GO_FILE', 'name', 'obo'])
    kws_set = set()

    def __init__(self):
        self.objdoc = DocOptParse(__doc__, self.kws_dict, self.kws_set)
        self.objsub = WrSubObo()

    def cli(self, prt=sys.stdout):
        """Command-line interface to print specified GO Terms from the DAG source ."""
        kws = self.objdoc.get_docargs(prt=None)
        # print("KWS", kws)
        goids = GetGOs().get_goids(kws.get('GO'), kws.get('GO_FILE'), sys.stdout)
        if not goids and 'name' in kws:
            goids = self.objsub.get_goids(kws['obo'], kws['name'])
        self.objsub.prt_goterms(kws['obo'], goids, prt, b_prt=False)
        print("Printing {N:6} GO IDs: {GOs}".format(N=len(goids), GOs=goids))


# Copyright (C) 2016-2018, DV Klopfenstein, H Tang. All rights reserved.
