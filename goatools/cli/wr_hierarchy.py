"""Command-line script to print a GO term's lower-level hierarchy.

Usage:
  wr_hier.py [GO ...] [options]

Options:
  -h --help      show this help message and exit

  -i <gofile.txt>  Read a file name containing a list of GO IDs
  --o=<outfile>    Output file in ASCII text format
  -f               Writes results to an ASCII file named after the GO term. e.g. hier_GO0002376.txt
  --up             Write report from GO term up to root

  --dag=<dag_file>    Ontologies in obo file [default: go-basic.obo].

  --gaf=<file.gaf>    Annotations from a gaf file
  --gene2go=<gene2go> Annotations from a gene2go file downloaded from NCBI

  --no_indent         Do not indent GO terms
  --max_indent=<int>  max indent depth for printing relative to GO Term
  --num_child=<int>   Print count of total number of children for each GO
  --short             If a branch has already been printed, do not re-print.
                      Print '===' instead of dashes to note the point of compression
  -r --relationship   Load and use the 'relationship' field
"""

from __future__ import print_function

__copyright__ = "Copyright (C) 2016-2018, DV Klopfenstein, H Tang. All rights reserved."
__author__ = "DV Klopfenstein"

import sys
from goatools.base import get_godag
from goatools.associations import get_tcntobj
from goatools.godag.obo_optional_attributes import OboOptionalAttrs
from goatools.cli.docopt_parse import DocOptParse
from goatools.cli.gos_get import GetGOs
from goatools.gosubdag.gosubdag import GoSubDag
from goatools.gosubdag.rpt.write_hierarchy import WrHierGO


def cli():
    """Command-line script to print a GO term's lower-level hierarchy."""
    objcli = WrHierCli(sys.argv[1:])
    fouts_txt = objcli.get_fouts()
    if fouts_txt:
        for fout_txt in fouts_txt:
            objcli.wrtxt_hier(fout_txt)
    else:
        objcli.prt_hier(sys.stdout)
    print(objcli.kws)

class WrHierCli(object):
    """Write hierarchy cli."""

    kws_set_all = set(['relationship', 'up', 'f'])
    kws_dct_all = set(['GO', 'dag', 'i', 'o', 'max_indent', 'num_child', 'no_indent', 'short',
                       'gaf', 'gene2go'])
    kws_dct_wr = set(['max_indent', 'num_child', 'no_indent', 'short', 'relationship'])

    def __init__(self, args=None, prt=sys.stdout):
        self.kws = DocOptParse(__doc__, self.kws_dct_all, self.kws_set_all).get_docargs(
            args, intvals=set(['max_indent', 'num_child']))
        opt_attrs = OboOptionalAttrs.attributes.intersection(self.kws.keys())
        godag = get_godag(self.kws['dag'], prt, optional_attrs=opt_attrs)
        # goids_usr = None if 'GO' not in self.kws else self.kws['GO']
        self.gosubdag = GoSubDag(godag.keys(), godag,
                                 relationships='relationship' in opt_attrs,
                                 tcntobj=get_tcntobj(godag, **self.kws),
                                 children=True,
                                 prt=prt)
        self.goids = GetGOs().get_goids(self.kws.get('GO'), self.kws.get('i'), sys.stdout)

    def get_fouts(self):
        """Get output filename."""
        fouts_txt = []
        if 'o' in self.kws:
            fouts_txt.append(self.kws['o'])
        if 'f' in self.kws:
            fouts_txt.append(self._get_fout_go())
        return fouts_txt

    def _get_fout_go(self):
        """Get the name of an output file based on the top GO term."""
        assert self.goids, "NO VALID GO IDs WERE PROVIDED"
        base = next(iter(self.goids)).replace(':', '')
        upstr = '_up' if 'up' in self.kws else ''
        return "hier_{BASE}{UP}.{EXT}".format(BASE=base, UP=upstr, EXT='txt')

    def wrtxt_hier(self, fout_txt):
        """Write hierarchy below specfied GO IDs to an ASCII file."""
        with open(fout_txt, 'wb') as prt:
            self.prt_hier(prt)
            print("  WROTE: {TXT}".format(TXT=fout_txt))

    def prt_hier(self, prt=sys.stdout):
        """Write hierarchy below specfied GO IDs."""
        objwr = WrHierGO(self.gosubdag, **self.kws)
        assert self.goids, "NO VALID GO IDs WERE PROVIDED"
        # kws = {k:v for k, v in self.kws.items() if k in self.kws_dct_wr}
        # objwr.write_hier_all(prt=prt, **kws)
                       # max_indent=None, num_child=None, short_prt=False):
        if 'up' not in objwr.usrset:
            for goid in self.goids:
                objwr.prt_hier_down(goid, prt)
        else:
            objwr.prt_hier_up(self.goids, prt)


# Copyright (C) 2016-2018, DV Klopfenstein, H Tang. All rights reserved.
