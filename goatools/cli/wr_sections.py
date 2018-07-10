"""Command-line interface to create an initial Python sections file

Usage:
  wr_sections.py [GO_FILE]
  wr_sections.py [GO_FILE] [options]

Options:
  -h --help                                 show this help message and exit

  -i <file.txt>, --ifile=<sections_in.txt>  Read or Write file name [default: sections_in.txt]
  -o <file.txt>, --ofile=<sections.txt>     write file name [default: sections.txt]
  --txt=<file.txt>                          Write file name [default: grouped_gos.txt]

  --py=<file.py>       Write the sections list into a Python file
  --xlsx=<file.xlsx>   Group user GO IDs and write the results into an xlsx file

  --obo=<file.obo>     Ontologies in obo file [default: go-basic.obo].
  --slims=<file.obo>   GO slims in obo file [default: goslim_generic.obo].

  --gaf=<file.gaf>     Annotations from a gaf file
  --gene2go=<gene2go>  Annotations from a gene2go file downloaded from NCBI

"""

from __future__ import print_function

__copyright__ = "Copyright (C) 2016-2018, DV Klopfenstein, H Tang. All rights reserved."
__author__ = "DV Klopfenstein"


import os
import sys

from goatools.base import get_godag
from goatools.associations import get_tcntobj

from goatools.cli.docopt_parse import DocOptParse
from goatools.cli.gos_get import GetGOs

from goatools.gosubdag.gosubdag import GoSubDag

from goatools.grouper.read_goids import read_sections
from goatools.grouper.grprdflts import GrouperDflts
from goatools.grouper.hdrgos import HdrgosSections
from goatools.grouper.grprobj import Grouper
from goatools.grouper.wr_sections import WrSectionsTxt
from goatools.grouper.wr_sections import WrSectionsPy
from goatools.grouper.sorter import Sorter
from goatools.grouper.wrxlsx import WrXlsxSortedGos


# pylint: disable=too-few-public-methods
class WrSectionsCli(object):
    """Class for command-line interface for creating GO term diagrams"""

    kws_dict = set(['GO_FILE', 'obo', 'slims',
                    'ifile', 'ofile', 'txt',
                    'py', 'xlsx',
                    'gaf', 'gene2go', 'taxid'])
    kws_set = set()

    def __init__(self, gosubdag=None):
        self.objdoc = DocOptParse(__doc__, self.kws_dict, self.kws_set)
        self.gosubdag = None if gosubdag is None else gosubdag

    def cli(self, prt=sys.stdout):
        """Command-line interface for go_draw script."""
        kws = self.objdoc.get_docargs(prt=None)
        godag = get_godag(kws['obo'], prt=None, loading_bar=False, optional_attrs=['relationship'])
        usrgos = GetGOs(godag, max_gos=200).get_usrgos(kws.get('GO_FILE'), prt)
        tcntobj = self._get_tcntobj(usrgos, godag, **kws)  # Gets TermCounts or None
        self.gosubdag = GoSubDag(usrgos, godag, relationships=True, tcntobj=tcntobj, prt=None)
        grprdflt = GrouperDflts(self.gosubdag, kws['slims'])
        ver_list = [godag.version, grprdflt.ver_goslims]
        prt.write("{VER}\n".format(VER="\n".join(ver_list)))
        sections = self._read_sections(kws['ifile'])
        # print("SECSECSEC", sections)
        hdrobj = HdrgosSections(self.gosubdag, grprdflt.hdrgos_dflt, sections)
        grprobj = Grouper("init", usrgos, hdrobj, self.gosubdag)
        # Write sections
        objsecwr = WrSectionsTxt(grprobj, ver_list)
        if not os.path.exists(kws['ifile']):
            objsecwr.wr_txt_section_hdrgos(kws['ifile'])
        objsecwr.wr_txt_section_hdrgos(kws['ofile'])
        objsecpy = WrSectionsPy(grprobj, ver_list)
        if 'py' in kws:
            objsecpy.wr_py_sections(kws['py'], sections, doc=godag.version)
        # Write user GO IDs in sections
        sortobj = Sorter(grprobj)
        objgowr = WrXlsxSortedGos("init", sortobj, ver_list)
        objgowr.wr_txt_gos(kws['txt'], sortby=objsecpy.fncsortnt)
        #objwr.wr_txt_section_hdrgos(kws['ofile'], sortby=objwr.fncsortnt)
        self._prt_cnt_usrgos(usrgos, sys.stdout)

    @staticmethod
    def _read_sections(ifile):
        """Read sections_in.txt file, if it exists."""
        if os.path.exists(ifile):
            return read_sections(ifile, exclude_ungrouped=True, prt=None)

    def _prt_cnt_usrgos(self, usrgos_read, prt):
        num_usrgos = len(self.gosubdag.go_sources)
        prt.write("{GOs:6} user GO IDs".format(GOs=num_usrgos))
        if len(usrgos_read) != num_usrgos:
            prt.write(" of {M} GO IDs read".format(M=len(usrgos_read)))
        prt.write("\n")

    @staticmethod
    def _get_tcntobj(goids, go2obj, **kws):
        """Get a TermCounts object if the user provides an annotation file, otherwise None."""
        # kws: gaf (gene2go taxid)
        if 'gaf' in kws or 'gene2go' in kws:
            # Get a reduced go2obj set for TermCounts
            _gosubdag = GoSubDag(goids, go2obj, rcntobj=False, prt=None)
            return get_tcntobj(_gosubdag.go2obj, **kws)  # TermCounts


# Copyright (C) 2016-2018, DV Klopfenstein, H Tang. All rights reserved.
