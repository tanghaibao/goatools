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

__copyright__ = (
    "Copyright (C) 2016-present, DV Klopfenstein, H Tang. All rights reserved."
)
__author__ = "DV Klopfenstein"


import os
import sys

from collections import namedtuple

from ..associations import get_tcntobj
from ..base import get_godag, logger
from ..godag.relationship_str import RelationshipStr
from ..gosubdag.gosubdag import GoSubDag
from ..gosubdag.rpt.wr_xlsx import GoDepth1LettersWr
from ..grouper.sorter import Sorter
from ..grouper.wrxlsx import WrXlsxSortedGos

from .docopt_parse import DocOptParse
from .gos_get import GetGOs
from .grouped import Grouped


# pylint: disable=too-few-public-methods
class CompareGOsCli:
    """Class for command-line interface for creating GO term diagrams"""

    kws_dict = set(
        [
            "GO_FILE",
            "sections",
            "S",
            "obo",
            "slims",
            "ofile",
            "xlsx",
            "gaf",
            "gene2go",
            "taxid",
        ]
    )
    kws_set = set(["verbose"])

    # Print fields to exclude, unless verbose is used
    excl_flds = {
        "level",
        "reldepth",
        "alt",
        "D1",
        "childcnt",
        "format_txt",
        "num_usrgos",
        "is_hdrgo",
        "is_usrgo",
        "hdr_idx",
        "hdr1usr01",
        "REL",
        "REL_short",
        "rel",
        "id",
    }

    def __init__(self, **kws):
        _objdoc = DocOptParse(__doc__, self.kws_dict, self.kws_set)
        self.kws = _objdoc.get_docargs(prt=None) if not kws else kws
        self.godag = get_godag(
            self.kws.get("obo"), prt=sys.stdout, optional_attrs=["relationship"]
        )
        _ini = _Init(self.godag)
        self.go_ntsets = _ini.get_go_ntsets(self.kws.get("GO_FILE"))
        self.go_all = set.union(*[nt.go_set for nt in self.go_ntsets])
        _tcntobj = _ini.get_tcntobj(self.go_all, **self.kws)  # Gets TermCounts or None
        self.gosubdag = GoSubDag(
            self.go_all, self.godag, True, tcntobj=_tcntobj, prt=sys.stdout
        )
        self.objgrpd = _ini.get_grouped(
            self.go_ntsets, self.go_all, self.gosubdag, **self.kws
        )
        # KWS: sortby hdrgo_sortby section_sortby

    def write(self, fout_xlsx=None, fout_txt=None, verbose=False):
        """Command-line interface for go_draw script."""
        sortby = self._get_fncsortnt(self.objgrpd.grprobj.gosubdag.prt_attr["flds"])
        kws_sort = {"sortby" if verbose else "section_sortby": sortby}
        sortobj = Sorter(self.objgrpd.grprobj, **kws_sort)
        # KWS: hdrgo_prt=True section_prt=None top_n=None use_sections=True
        # RET: {sortobj, sections, hdrgo_prt} or {sortobj flat hdrgo_prt}
        desc2nts = sortobj.get_desc2nts_fnc(
            hdrgo_prt=verbose, section_prt=True, top_n=None, use_sections=True
        )
        # print('FFFF', desc2nts['flds'])
        # Write user GO IDs in sections
        objgowr = WrXlsxSortedGos("init", sortobj, self.objgrpd.ver_list)
        if fout_xlsx is not None:
            kws_xlsx = {"shade_hdrgos": verbose}
            if not verbose:
                kws_xlsx["prt_flds"] = [
                    f for f in desc2nts["flds"] if f not in self.excl_flds
                ]
            self._adj_hdrs(kws_xlsx, desc2nts)
            objgowr.wr_xlsx_nts(fout_xlsx, desc2nts, **kws_xlsx)
            fout_desc = "{BASE}_desc.txt".format(BASE=os.path.splitext(fout_xlsx)[0])
            self._wr_ver_n_key(fout_desc, verbose)
        if fout_txt is not None:
            self._wr_txt_nts(fout_txt, desc2nts, objgowr, verbose)
        if fout_xlsx is None and fout_txt is None:
            self._prt_ver_n_key(sys.stdout, verbose)
            prtfmt = self._get_prtfmt(objgowr, verbose)
            summary_dct = objgowr.prt_txt_desc2nts(sys.stdout, desc2nts, prtfmt)
            self._prt_ver_n_key(sys.stdout, verbose)
            if summary_dct:
                print(
                    "\n{N} GO IDs in {S} sections".format(
                        N=desc2nts["num_items"], S=desc2nts["num_sections"]
                    )
                )

    def _adj_hdrs(self, kws_xlsx, desc2nts):
        """Replace xlsx column header, fileN, with base input filenames"""
        filehdrs = [nt.hdr for nt in self.go_ntsets]
        num_files = len(filehdrs)
        if num_files == len(set(filehdrs)):
            kws_xlsx["hdrs"] = filehdrs + list(desc2nts["flds"][num_files:])

    def _get_prtfmt(self, objgowr, verbose):
        """Get print format containing markers."""
        prtfmt = objgowr.get_prtfmt("fmt")
        prtfmt = prtfmt.replace("# ", "")
        if not verbose:
            prtfmt = prtfmt.replace("{hdr1usr01:2}", "")
            prtfmt = prtfmt.replace("{childcnt:3} L{level:02} ", "")
            prtfmt = prtfmt.replace("{num_usrgos:>4} uGOs ", "")
            prtfmt = prtfmt.replace("{D1:5} {REL} {rel}", "")
            prtfmt = prtfmt.replace("R{reldepth:02} ", "")
        marks = "".join(["{{{}}}".format(nt.fileN) for nt in self.go_ntsets])
        return "{MARKS} {PRTFMT}".format(MARKS=marks, PRTFMT=prtfmt)

    @staticmethod
    def _get_fncsortnt(flds):
        """Return a sort function for sorting header GO IDs found in sections."""
        if "tinfo" in flds:
            return lambda ntgo: [ntgo.NS, -1 * ntgo.tinfo, ntgo.depth, ntgo.alt]
        if "dcnt" in flds:
            return lambda ntgo: [ntgo.NS, -1 * ntgo.dcnt, ntgo.depth, ntgo.alt]
        return lambda ntgo: [ntgo.NS, -1 * ntgo.depth, ntgo.alt]

    def _wr_txt_nts(self, fout_txt, desc2nts, objgowr, verbose):
        """Write grouped and sorted GO IDs to GOs."""
        with open(fout_txt, "w", encoding="utf-8") as prt:
            self._prt_ver_n_key(prt, verbose)
            prt.write("\n\n")
            prt.write(
                "# ----------------------------------------------------------------\n"
            )
            prt.write("# - Sections and GO IDs\n")
            prt.write(
                "# ----------------------------------------------------------------\n"
            )
            prtfmt = self._get_prtfmt(objgowr, verbose)
            summary_dct = objgowr.prt_txt_desc2nts(prt, desc2nts, prtfmt)
            if summary_dct:
                print(
                    "  {N:>5} GO IDs WROTE: {FOUT} ({S} sections)".format(
                        N=desc2nts["num_items"],
                        FOUT=fout_txt,
                        S=desc2nts["num_sections"],
                    )
                )
            else:
                print("  WROTE: {TXT}".format(TXT=fout_txt))

    def _wr_ver_n_key(self, fout_txt, verbose):
        """Write GO DAG version and key indicating presence of GO ID in a list."""
        with open(fout_txt, "w", encoding="utf-8") as prt:
            self._prt_ver_n_key(prt, verbose)
            print("               WROTE: {TXT}".format(TXT=fout_txt))

    def _prt_ver_n_key(self, prt, verbose):
        """Print GO DAG version and key indicating presence of GO ID in a list."""
        pre = "# "
        prt.write(
            "# ----------------------------------------------------------------\n"
        )
        prt.write("# - Description of GO ID fields\n")
        prt.write(
            "# ----------------------------------------------------------------\n"
        )
        prt.write(
            "# Versions:\n#    {VER}\n".format(
                VER="\n#    ".join(self.objgrpd.ver_list)
            )
        )
        prt.write("\n# Marker keys:\n")
        for ntgos in self.go_ntsets:
            prt.write("#     X -> GO is present in {HDR}\n".format(HDR=ntgos.hdr))
        if verbose:
            prt.write("\n# Markers for header GO IDs and user GO IDs:\n")
            prt.write("#     '**' -> GO term is both a header and a user GO ID\n")
            prt.write("#     '* ' -> GO term is a header, but not a user GO ID\n")
            prt.write("#     '  ' -> GO term is a user GO ID\n")
        prt.write("\n# GO Namspaces:\n")
        prt.write("#     BP -> Biological Process\n")
        prt.write("#     MF -> Molecular Function\n")
        prt.write("#     CC -> Cellular Component\n")
        if verbose:
            prt.write("\n# Example fields: 5 uGOs   362  47 L04 D04 R04\n")
            prt.write(
                "#     N uGOs         -> number of user GO IDs under this GO header\n"
            )
            prt.write("#     First integer  -> number of GO descendants\n")
            prt.write(
                "#     Second integer -> number of GO children for the current GO ID\n"
            )
        prt.write("\n# Depth information:\n")
        if not verbose:
            prt.write("#     int -> number of GO descendants\n")
        if verbose:
            prt.write("#     Lnn -> level (minimum distance from root to node)\n")
        prt.write("#     Dnn -> depth (maximum distance from root to node)\n")
        if verbose:
            prt.write("#     Rnn -> depth accounting for relationships\n\n")
            RelationshipStr().prt_keys(prt, pre)
        if verbose:
            prt.write("\n")
            objd1 = GoDepth1LettersWr(self.gosubdag.rcntobj)
            objd1.prt_header(prt, "DEPTH-01 GO terms and their aliases", pre)
            objd1.prt_txt(prt, pre)


class _Init:
    """Initialize object."""

    def __init__(self, godag):
        self.godag = godag

    def get_tcntobj(self, go_all, **kws):
        """Get a TermCounts object if the user provides an annotation file, otherwise None."""
        # kws: gaf (gene2go taxid)
        if "gaf" in kws or "gene2go" in kws:
            # Get a reduced go2obj set for TermCounts
            _gosubdag = GoSubDag(go_all, self.godag, rcntobj=False, prt=None)
            return get_tcntobj(_gosubdag.go2obj, **kws)  # TermCounts
        return None

    def get_grouped(self, go_ntsets, go_all, gosubdag, **kws):
        """Get Grouped object."""
        kws_grpd = {k: v for k, v in kws.items() if k in Grouped.kws_dict}
        kws_grpd["go2nt"] = self._init_go2ntpresent(go_ntsets, go_all, gosubdag)
        return Grouped(gosubdag, self.godag.version, **kws_grpd)

    @staticmethod
    def _init_go2ntpresent(go_ntsets, go_all, gosubdag):
        """Mark all GO IDs with an X if present in the user GO list."""
        go2ntpresent = {}
        flds = " ".join(nt.fileN for nt in go_ntsets)
        ntobj = namedtuple("NtPresent", flds)
        # Get present marks for GO sources
        for goid_all in go_all:
            present_true = [goid_all in nt.go_set for nt in go_ntsets]
            present_str = ["X" if tf else "." for tf in present_true]
            go2ntpresent[goid_all] = ntobj._make(present_str)
        # Get present marks for all other GO ancestors
        goids_ancestors = set(gosubdag.go2obj).difference(go2ntpresent)
        assert not goids_ancestors.intersection(go_all)
        strmark = ["." for _ in range(len(go_ntsets))]
        for goid in goids_ancestors:
            go2ntpresent[goid] = ntobj._make(strmark)
        return go2ntpresent

    def get_go_ntsets(self, go_fins):
        """For each file containing GOs, extract GO IDs, store filename and header."""
        nts = []
        go_fins = list(go_fins)
        ntobj = namedtuple("NtGOFiles", "fileN hdr go_set go_fin")
        go_sets = self._init_go_sets(go_fins)
        hdrs = [os.path.splitext(os.path.basename(f))[0] for f in go_fins]
        assert len(go_fins) == len(go_sets)
        assert len(go_fins) == len(hdrs)
        goids = set()
        for idx, (hdr, go_set, go_fin) in enumerate(zip(hdrs, go_sets, go_fins), 1):
            goids.update(go_set)
            if not go_set:
                logger.warning("NO GO IDs FOUND IN %s", go_fin)
            nts.append(
                ntobj(
                    fileN="file{I}".format(I=idx), hdr=hdr, go_set=go_set, go_fin=go_fin
                )
            )
        if not goids:
            logger.warning("NO GO IDs FOUND")
            sys.exit(1)
        return nts

    def _init_go_sets(self, go_fins):
        """Get lists of GO IDs."""
        go_sets = []
        assert go_fins, "EXPECTED FILES CONTAINING GO IDs"
        assert len(go_fins) >= 2, "EXPECTED 2+ GO LISTS. FOUND: {L}".format(
            L=" ".join(go_fins)
        )
        obj = GetGOs(self.godag)
        for fin in go_fins:
            assert os.path.exists(fin), "GO FILE({F}) DOES NOT EXIST".format(F=fin)
            go_sets.append(obj.get_usrgos(fin, sys.stdout))
        return go_sets


# Copyright (C) 2016-present, DV Klopfenstein, H Tang. All rights reserved.
