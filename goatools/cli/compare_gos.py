"""Compare two or more sets of GO IDs. Best done using sections.

Usage:
  compare_gos.py [GO_FILE] ...
  compare_gos.py [GO_FILE] ... [options]

Options:
  -h --help            show this help message and exit

  -s <sections.txt> --sections=<sections.txt>  Sections file for grouping
  -S <sections module str>                     Python module with SECTIONS variable

  -o <file.txt>, --ofile=<file.txt>    write comparison of GO IDs into ASCII file
  --xlsx=<file.xlsx>   write comparison of GO IDs into an xlsx file
  -v --verbose         Print sections as GO headers followed by each header's user GOs

  --obo=<file.obo>     Ontologies in obo file [default: go-basic.obo].
  --slims=<file.obo>   GO slims in obo file [default: goslim_generic.obo].

  --gaf=<file.gaf>     Annotations from a gaf file
  --gene2go=<gene2go>  Annotations from a gene2go file downloaded from NCBI

"""

from __future__ import print_function

__copyright__ = "Copyright (C) 2016-2019, DV Klopfenstein, H Tang. All rights reserved."
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
from goatools.grouper.sorter import Sorter
from goatools.grouper.wrxlsx import WrXlsxSortedGos


# pylint: disable=too-few-public-methods
class CompareGOsCli(object):
    """Class for command-line interface for creating GO term diagrams"""

    kws_dict = set(['GO_FILE',
                    'sections', 'S',
                    'obo', 'slims',
                    'ofile', 'xlsx',
                    'gaf', 'gene2go', 'taxid',
                   ])
    kws_set = set(['verbose'])

    # Print fields to exclude, unless verbose is used
    excl_flds = {'level', 'reldepth', 'alt', 'D1', 'childcnt',
                 'format_txt', 'num_usrgos', 'is_hdrgo', 'is_usrgo', 'hdr_idx', 'hdr1usr01',
                 'REL', 'REL_short', 'rel', 'id'}

    def __init__(self, **kws):
        _objdoc = DocOptParse(__doc__, self.kws_dict, self.kws_set)
        self.kws = _objdoc.get_docargs(prt=None) if not kws else kws
        # for key, val in self.kws.items():
        #     print('WWWWWWWWWWWWWWWWWWW', key, val)
        self.godag = get_godag(self.kws.get('obo'), prt=sys.stdout,
                               loading_bar=False, optional_attrs=['relationship'])
        self.go_fins = self.kws.get('GO_FILE')
        self.go_sets = self._init_go_sets()
        self.go_all = set.union(*self.go_sets)
        self.objgrpd = self._init_grouped()
        # KWS: sortby hdrgo_sortby section_sortby

    def write(self, fout_xlsx=None, fout_txt=None, verbose=False):
        """Command-line interface for go_draw script."""
        #print('VVVVVVVVVVVVVVVV verbose', verbose)
        sys.stdout.write("{VER}\n".format(VER="\n".join(self.objgrpd.ver_list)))
        sortby = self._get_fncsortnt(self.objgrpd.grprobj.gosubdag.prt_attr['flds'])
        kws_sort = {'sortby' if verbose else 'section_sortby': sortby}
        sortobj = Sorter(self.objgrpd.grprobj, **kws_sort)
        # KWS: hdrgo_prt=True section_prt=None top_n=None use_sections=True
        # RET: {sortobj, sections, hdrgo_prt} or {sortobj flat hdrgo_prt}
        desc2nts = sortobj.get_desc2nts_fnc(
            hdrgo_prt=verbose,
            section_prt=True,
            top_n=None,
            use_sections=True)
        print('FFFF', desc2nts['flds'])
        # Write user GO IDs in sections
        objgowr = WrXlsxSortedGos("init", sortobj, self.objgrpd.ver_list)
        if fout_xlsx is not None:
            kws_xlsx = {'shade_hdrgos':verbose}
            if not verbose:
                kws_xlsx['prt_flds'] = [f for f in desc2nts['flds'] if f not in self.excl_flds]
            objgowr.wr_xlsx_nts(fout_xlsx, desc2nts, **kws_xlsx)
        if fout_txt is not None:
            prtfmt = objgowr.get_prtfmt('fmt')
            print('FFFFFFFFFFFF', prtfmt)
            objgowr.wr_txt_nts(fout_txt, desc2nts, prtfmt=prtfmt)
        if fout_xlsx is None and fout_txt is None:
            summary_dct = objgowr.prt_txt_desc2nts(sys.stdout, desc2nts, prtfmt=None)
            if summary_dct:
                print(sortobj.grprobj.fmtsum.format(ACTION='PRINTED:', FILE='', **summary_dct))
        # SUMMARY: hdr GOs(24 in 15 sections, N/A unused) READ: data/compare_gos/sections.txt
        self._prt_cnt_usrgos(self.go_all, sys.stdout)

    @staticmethod
    def _get_fncsortnt(flds):
        """Return a sort function for sorting header GO IDs found in sections."""
        if 'tinfo' in flds:
            return lambda ntgo: [ntgo.NS, -1*ntgo.tinfo, ntgo.depth, ntgo.alt]
        if 'dcnt' in flds:
            return lambda ntgo: [ntgo.NS, -1*ntgo.dcnt, ntgo.depth, ntgo.alt]
        return lambda ntgo: [ntgo.NS, -1*ntgo.depth, ntgo.alt]

    def _init_go_sets(self, prt=sys.stdout):
        """Get lists of GO IDs."""
        go_sets = []
        assert self.go_fins, "EXPECTED FILES CONTAINING GO IDs"
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
        hdrs = [os.path.splitext(os.path.basename(f))[0] for f in self.go_fins]
        ntobj = namedtuple('NtPresent', " ".join(hdrs))
        for goid_all in self.go_all:
            present_true = [goid_all in gos for gos in self.go_sets]
            present_str = ['X' if tf else '.' for tf in present_true]
            go2present[goid_all] = ntobj._make(present_str)
        return go2present

# Copyright (C) 2016-2019, DV Klopfenstein, H Tang. All rights reserved.
