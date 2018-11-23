"""Compare two or more sets of GO IDs. Best done using sections.

Usage:
  compare_gos.py [GO_FILE] ...
  compare_gos.py [GO_FILE] ... [options]

Options:
  -h --help                                 show this help message and exit

  -s <sections.txt> --sections=<sections.txt>  Sections file for grouping
  -S <sections module str>             Sections file for grouping

  -o <file.txt>, --ofile=<file.txt>    write comparison of GO IDs into ASCII file
  --xlsx=<file.xlsx>   write comparison of GO IDs into an xlsx file

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
from collections import namedtuple
# from collections import OrderedDict

from goatools.base import get_godag
from goatools.associations import get_tcntobj

from goatools.cli.docopt_parse import DocOptParse
from goatools.cli.gos_get import GetGOs

from goatools.gosubdag.gosubdag import GoSubDag

from goatools.cli.grouped import Grouped
#### from goatools.grouper.read_goids import read_sections
#### from goatools.grouper.grprdflts import GrouperDflts
#### from goatools.grouper.hdrgos import HdrgosSections
#### from goatools.grouper.grprobj import Grouper
#### from goatools.grouper.wr_sections import WrSectionsTxt
#### from goatools.grouper.wr_sections import WrSectionsPy
from goatools.grouper.wr_sections import get_fncsortnt
from goatools.grouper.sorter import Sorter
from goatools.grouper.wrxlsx import WrXlsxSortedGos


# pylint: disable=too-few-public-methods
class CompareGOs(object):
    """Class for command-line interface for creating GO term diagrams"""

    kws_dict = set(['GO_FILE',
                    'sections', 'S',
                    'obo', 'slims',
                    'ofile', 'xlsx',
                    'gaf', 'gene2go', 'taxid'])
    kws_set = set()

    def __init__(self, **kws):
        _objdoc = DocOptParse(__doc__, self.kws_dict, self.kws_set)
        self.kws = _objdoc.get_docargs(prt=None) if not kws else kws
        self.godag = get_godag(self.kws.get('obo'), prt=sys.stdout,
                               loading_bar=False, optional_attrs=['relationship'])
        self.go_fins = self.kws.get('GO_FILE')
        self.go_sets = self._init_go_sets()
        self.go_all = set.union(*self.go_sets)
        self.objgrpd = self._init_grouped()
        self.sortobj_all = Sorter(self.objgrpd.grprobj)

    def cli(self, prt=sys.stdout):
        """Command-line interface for go_draw script."""
#        self.gosubdag = GoSubDag(usrgos, godag, relationships=True, tcntobj=tcntobj, prt=None)
        prt.write("{VER}\n".format(VER="\n".join(self.objgrpd.ver_list)))
#        # print("SECSECSEC", sections)
#        grprobj = Grouper("init", usrgos, hdrobj, self.gosubdag)
        ### Write sections
        ## objsecwr = WrSectionsTxt(self.grprobj, ver_list)
        ## if not os.path.exists(kws['ifile']):
        ##     objsecwr.wr_txt_section_hdrgos(kws['ifile'])
        ## objsecwr.wr_txt_section_hdrgos(kws['ofile'])
        ## objsecpy = WrSectionsPy(self.grprobj, self.ver_list)
        ## if 'py' in kws:
        ##     objsecpy.wr_py_sections(kws['py'], sections, doc=godag.version)
        # Write user GO IDs in sections
        objgowr = WrXlsxSortedGos("init", self.sortobj_all, self.objgrpd.ver_list)
        sortby = get_fncsortnt(self.objgrpd.grprobj.gosubdag.prt_attr['flds'])
        objgowr.wr_txt_gos(self.kws['ofile'], sortby=sortby)
        if 'xlsx' in self.kws:
            objgowr.wr_xlsx_gos(self.kws['xlsx'], sortby=sortby)
        #objwr.wr_txt_section_hdrgos(kws['ofile'], sortby=objwr.fncsortnt)
        # hdr GOs(24 in 15 sections, N/A unused) READ:  data/compare_gos/sections.txt
        self._prt_cnt_usrgos(self.go_all, sys.stdout)

    def _init_go_sets(self, prt=sys.stdout):
        """Get lists of GO IDs."""
        go_sets = []
        assert len(self.go_fins) >= 2, "EXPECTED 2+ GO LISTS. FOUND: {L}".format(
            L=' '.join(self.go_fins))
        obj = GetGOs(self.godag)
        for fin in self.go_fins:
            assert os.path.exists(fin), "GO FILE({F}) DOES NOT EXIST".format(F=fin)
            go_sets.append(obj.get_usrgos(fin, prt))
        return go_sets

    def _prt_cnt_usrgos(self, usrgos_read, prt):
        num_usrgos = len(self.objgrpd.gosubdag.go_sources)
        prt.write("{GOs:6} user GO IDs".format(GOs=num_usrgos))
        if len(usrgos_read) != num_usrgos:
            prt.write(" of {M} GO IDs read".format(M=len(usrgos_read)))
        prt.write("\n")

    def _get_tcntobj(self, **kws):
        """Get a TermCounts object if the user provides an annotation file, otherwise None."""
        # kws: gaf (gene2go taxid)
        if 'gaf' in kws or 'gene2go' in kws:
            # Get a reduced go2obj set for TermCounts
            _gosubdag = GoSubDag(self.go_all, self.godag, rcntobj=False, prt=None)
            return get_tcntobj(_gosubdag.go2obj, **kws)  # TermCounts

    def _init_grouped(self):
        """Get Grouped object."""
        _tcntobj = self._get_tcntobj(**self.kws)  # Gets TermCounts or None
        kws_grpd = {k:v for k, v in self.kws.items() if k in Grouped.kws_dict}
        kws_grpd['go2nt'] = self._init_go2present()
        return Grouped(self.go_all, self.godag, _tcntobj, **kws_grpd)

    def _init_go2present(self):
        """Mark all GO IDs set membership in each section."""
        go2present = {}
        hdrs = ["i{N}".format(N=n) for n in range(len(self.go_sets))]
        ntobj = namedtuple('NtPresent', " ".join(hdrs))
        for goid_all in self.go_all:
            go2present[goid_all] = ntobj._make([goid_all in gos for gos in self.go_sets])
        return go2present

# Copyright (C) 2016-2018, DV Klopfenstein, H Tang. All rights reserved.
