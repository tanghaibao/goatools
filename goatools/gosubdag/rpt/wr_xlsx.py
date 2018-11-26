"""Used to operate on a sub-graph of a larger GO DAG."""

__copyright__ = "Copyright (C) 2016-2019, DV Klopfenstein, H Tang, All rights reserved."
__author__ = "DV Klopfenstein"

import sys
import re
import collections as cx
from goatools.gosubdag.gosubdag import GoSubDag
from goatools.wr_tbl import wr_xlsx as wr_xlsx_tbl
from goatools.wr_tbl import wr_xlsx_sections as wr_xlsx_sections_tbl
from goatools.wr_tbl import get_lines
#### from goatools.wr_tbl import prt_txt

class GoSubDagWr(object):
    """Contains a sub-graph of the original obo from geneontology.org."""

    fld2col_widths = {
        'NS' : 3,
        'dcnt' : 6,
        'level' : 4,
        'depth' : 4,
        'GO' : 12,
        'D1' : 6,
        'GO_name' : 45}

    def __init__(self, go2obj):
        self.go2obj = go2obj

    def wr_xlsx(self, fout_xlsx, goids, sortby=None, **kws_usr):
        """Write goids into a table."""
        nts = GoSubDag(goids, self.go2obj).get_nts(goids, sortby)
        kws_wr = kws_usr.copy()
        if 'fld2col_widths' not in kws_wr:
            kws_wr['fld2col_widths'] = self.fld2col_widths
        wr_xlsx_tbl(fout_xlsx, nts, **kws_wr)

    def wr_xlsx_sections(self, fout_xlsx, sections, sortby=None, **kws_usr):
        """Write goids into a table."""
        nts = self.get_nts_sections(sections, sortby)
        kws_wr = kws_usr.copy()
        if 'fld2col_widths' not in kws_wr:
            kws_wr['fld2col_widths'] = self.fld2col_widths
        else:
            fld2col_widths = self.fld2col_widths.copy()
            for fld, wid in kws_usr['fld2col_widths'].items():
                fld2col_widths[fld] = wid
            kws_wr['fld2col_widths'] = fld2col_widths
        wr_xlsx_sections_tbl(fout_xlsx, nts, **kws_wr)

    def get_nts_sections(self, sections, sortby=None):
        """Given a list of sections containing GO IDs, get a list of sections w/GO nts."""
        goids = self.get_goids_sections(sections)
        gosubdag = GoSubDag(goids, self.go2obj)
        return [(sec, gosubdag.get_nts(gos, sortby)) for sec, gos in sections]

    @staticmethod
    def get_goids_sections(sections):
        """Return all the GO IDs in a 2-D sections list."""
        goids_all = set()
        for _, goids_sec in sections:
            goids_all |= set(goids_sec)
        return goids_all


def read_d1_letter(fin_txt):
    """Reads letter aliases from a text file created by GoDepth1LettersWr."""
    go2letter = {}
    re_goid = re.compile(r"(GO:\d{7})")
    with open(fin_txt) as ifstrm:
        for line in ifstrm:
            mtch = re_goid.search(line)
            if mtch and line[:1] != ' ':
                # Alias is expected to be the first character
                go2letter[mtch.group(1)] = line[:1]
    return go2letter

class GoDepth1LettersWr(object):
    """Writes reports for a GoDepth1Letters object."""

    str2ns = {'biological_process': 'BP', 'molecular_function': 'MF', 'cellular_component': 'CC'}
    hdrs = ['D1', 'NS', 'descendants', 'depth', 'GO', 'GO description']

    def __init__(self, rcntobj):
        self.ns2nt = self._init_ns2nt(rcntobj)
        self.goone2ntletter = rcntobj.goone2ntletter

    def prt_txt(self, prt=sys.stdout, pre=''):
        """Print letters, descendant count, and GO information."""
        data_nts = self.get_d1nts()
        for ntdata in data_nts:
            prt.write("{PRE}{L:1} {NS} {d:6,} D{D:02} {GO} {NAME}\n".format(
                PRE=pre,
                L=ntdata.D1,
                d=ntdata.dcnt,
                NS=ntdata.NS,
                D=ntdata.depth,
                GO=ntdata.GO,
                NAME=ntdata.name))
        return data_nts

    def wr_xlsx(self, fout_xlsx="gos_depth01.xlsx", **kws):
        """Write xlsx table of depth-01 GO terms and their letter representation."""
        data_nts = self.get_d1nts()
        if 'fld2col_widths' not in kws:
            kws['fld2col_widths'] = {'D1': 6, 'NS':3, 'depth': 5, 'GO': 12, 'name': 40}
        if 'hdrs' not in kws:
            kws['hdrs'] = self.hdrs
        wr_xlsx_tbl(fout_xlsx, data_nts, **kws)

    def wr_txt(self, fout_txt="gos_depth01.txt", title=None):
        """write text table of depth-01 GO terms and their letter representation."""
        with open(fout_txt, 'w') as prt:
            self.prt_header(prt, title)
            data_nts = self.prt_txt(prt)
            sys.stdout.write("  {N:>5} items WROTE: {TXT}\n".format(
                N=len(data_nts), TXT=fout_txt))

    @staticmethod
    def prt_header(prt, title=None, pre=''):
        """write text table of depth-01 GO terms and their letter representation."""
        if title is not None:
            prt.write("{PRE}{TITLE}\n".format(TITLE=title, PRE=pre))
            prt.write('{PRE}\n'.format(PRE=pre))
        prt.write("{PRE}    D1 : Letter representing the depth-01 GO term\n".format(PRE=pre))
        prt.write("{PRE}    dcnt: Total number of all descendants\n".format(PRE=pre))
        prt.write("{PRE}    dep: Depth; The maximum length path to ".format(PRE=pre))
        prt.write("{PRE}leaf-level (childless) GO descendant(s)\n".format(PRE=pre))
        prt.write("{PRE}\n".format(PRE=pre))
        prt.write("{PRE}D1 NS  dcnt dep GO ID      Description\n".format(PRE=pre))
        prt.write("{PRE}- -- ------ --- ---------- ------------------------------\n".format(
            PRE=pre))

    def wr_tex(self, fout_tex="gos_depth01.tex"):
        """write text table of depth-01 GO terms and their letter representation."""
        data_nts = self.get_d1nts()
        joinchr = " & "
        #pylint: disable=anomalous-backslash-in-string
        eol = " \\\\\n"
        with open(fout_tex, 'w') as prt:
            prt.write("\\begin{table}[!ht]\n")
            prt.write("\\begin{tabular}{|p{.5cm} | p{.5cm} | >{\\raggedleft}p{.9cm} ")
            prt.write("|p{.7cm} |p{1.8cm} |p{9cm}|}\n")
            prt.write("\multicolumn{6}{c}{} \\\\\n")
            prt.write("\hline\n")
            prt.write("\\rowcolor{gray!10}\n")
            prt.write("{HDRS}{EOL}".format(
                HDRS=joinchr.join(next(iter(data_nts))._fields), EOL=eol))
            prt.write("\hline\n")
            for idx, line in enumerate(get_lines(data_nts, joinchr=joinchr, eol=eol)):
                if idx%2 == 1:
                    prt.write("\\rowcolor{gray!7}\n")
                line.replace('_', '\\_')
                prt.write(line)
            prt.write("\hline\n")
            prt.write("\end{tabular}\n")
            caption = ("The descendant counts of GO terms at depth-01 are highly skewed. The "
                       "root term, \textit{biological\_process} has over twenty GO children at "
                       "depth-01 shown in the table sorted by their number of descendants "
                       "(dcnt) with \textit{cellular process} at the top having 18k+ "
                       "descendants and \textit{cell killing} near the bottom having only "
                       "about 100 descendants. The first column (D1) contains a letter used as "
                       "an alias for each depth-01 GO term. The second column represents the "
                       "number of descendants from the specified GO term from down to the total  "
                       "of its descendant leaf-level GO terms, which have no child GO terms.")
            prt.write("\caption{{{TEXT}}}\n\n".format(TEXT=caption))
            prt.write("\label{table:supptbl_d1}\n")
            prt.write("\end{table}\n")
            sys.stdout.write("  {N:>5} items WROTE: {TXT}\n".format(
                N=len(data_nts), TXT=fout_tex))

    def get_d1nts(self):
        """Get letters for depth-01 GO terms, descendants count, and GO information."""
        data = []
        ntdata = cx.namedtuple("NtPrt", "D1 NS dcnt depth GO name")
        namespace = None
        for ntlet in sorted(self.goone2ntletter.values(),
                            key=lambda nt: [nt.goobj.namespace, -1 * nt.dcnt, nt.D1]):
            goobj = ntlet.goobj
            goid = goobj.id
            assert len(goobj.parents) == 1
            if namespace != goobj.namespace:
                namespace = goobj.namespace
                ntns = self.ns2nt[namespace]
                pobj = ntns.goobj
                ns2 = self.str2ns[goobj.namespace]
                data.append(ntdata._make([" ", ns2, ntns.dcnt, pobj.depth, pobj.id, pobj.name]))
            data.append(ntdata._make(
                [ntlet.D1, self.str2ns[namespace], ntlet.dcnt, goobj.depth, goid, goobj.name]))
        return data

    @staticmethod
    def _init_ns2nt(rcntobj):
        """Save depth-00 GO terms ordered using descendants cnt."""
        go2dcnt = rcntobj.go2dcnt
        ntobj = cx.namedtuple("NtD1", "D1 dcnt goobj")
        d0s = rcntobj.depth2goobjs[0]
        ns_nt = [(o.namespace, ntobj(D1="", dcnt=go2dcnt[o.id], goobj=o)) for o in d0s]
        return cx.OrderedDict(ns_nt)

# Copyright (C) 2016-2019, DV Klopfenstein, H Tang, All rights reserved.
