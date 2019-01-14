"""Report the level/depth summaries of all GO terms or a subset of GO terms.

   Level is the length of a shortest path to a GO term.
   Depth is the length of a longest path to a GO term.

    Example:

    >>> godag = GODag("go-basic.obo")
    >>> reporter = RptLevDepth(godag)
    >>> reporter.write_summary_cnts_all()
        go-basic.obo: format-version(1.2) data-version(releases/2016-03-01)
        46162 nodes imported
        Dep <-Depth Counts->  <-Level Counts->
        Lev   BP    MF    CC    BP    MF    CC
        --- ----  ----  ----  ----  ----  ----
        00     1     1     1     1     1     1
        01    24    19    24    24    19    24
        02   126   131   195   222   154   336
        03   955   493   502  1913   739  1145
        04  1950  1464   566  4506  1812  1300
        05  3371  3765   994  6958  3985   763
        06  4275  1784   748  6945  1888   268
        07  4634  1005   549  4890   903    53
        08  4125   577   216  2022   352     6
        09  3478   315    86   739   109     1
        10  2364   164    12   188    38     0
        11  1583   170     3    37    21     0
        12  1014    68     1     1     0     0
        13   400    49     0     0     0     0
        14   107    13     0     0     0     0
        15    28     3     0     0     0     0
        16    11     0     0     0     0     0

"""

__copyright__ = "Copyright (C) 2015-2018, DV Klopfenstein, H Tang, All rights reserved."
__author__ = "DV Klopfenstein"

import sys
import collections as cx

def prt_lev_depth(godag, prt=sys.stdout):
    """Write the counts of GO terms found at each level and depth of the GO DAG."""
    RptLevDepth(godag, prt).write_summary_cnts_all()

class RptLevDepth(object):
    """Reports a summary of GO term levels in depths."""

    nss = ['biological_process', 'molecular_function', 'cellular_component']

    def __init__(self, obodag, log=sys.stdout):
        self.obo = obodag
        self.log = log
        self.title = "GO Counts in {VER}".format(VER=self.obo.version)

    def wr_xlsx(self, fout_xlsx):
        """Write counts of GO terms at all levels and depths."""
        from goatools.wr_tbl import wr_xlsx
        kws = {
            'title' : self.title,
            'hdrs' :["Dep/Lev",
                     "BP Depth", "MF Depth", "CC Depth",
                     "BP Level", "MF Level", "CC Level"]}
        wr_xlsx(fout_xlsx, self.get_data(), **kws)

    def wr_txt(self, fout_txt):
        """Write counts of GO terms at all levels and depths."""
        from goatools.wr_tbl import prt_txt
        data = self.get_data()
        with open(fout_txt, 'w') as prt:
            prtfmt = "{Depth_Level:>7} " \
                     "{BP_D:6,} {MF_D:6,} {CC_D:>6,} " \
                     "{BP_L:>6,} {MF_L:>6,} {CC_L:>6,}\n"
            prt.write("{TITLE}\n\n".format(TITLE=self.title))
            prt.write("         |<---- Depth ---->|  |<---- Level ---->|\n")
            prt.write("Dep/Lev     BP     MF     CC     BP     MF     CC\n")
            prt.write("-------  -----  -----  -----  -----  -----  -----\n")
            prt_txt(prt, data, prtfmt=prtfmt, title=self.title)
            sys.stdout.write("  {N:>5,} items WROTE: {TXT}\n".format(
                N=len(data), TXT=fout_txt))

    def prttex_summary_cnts_all(self, prt=sys.stdout):
        """Print LaTeX format summary of level and depth counts for all active GO Terms."""
        cnts = self.get_cnts_levels_depths_recs(set(self.obo.values()))
        self._prttex_summary_cnts(prt, cnts)

    def write_summary_cnts_all(self):
        """Write summary of level and depth counts for all active GO Terms."""
        cnts = self.get_cnts_levels_depths_recs(set(self.obo.values()))
        self._write_summary_cnts(cnts)

    def write_summary_cnts(self, go_ids):
        """Write summary of level and depth counts for specific GO ids."""
        obo = self.obo
        cnts = self.get_cnts_levels_depths_recs([obo.get(GO) for GO in go_ids])
        self._write_summary_cnts(cnts)

    def write_summary_cnts_goobjs(self, goobjs):
        """Write summary of level and depth counts for active GO Terms."""
        cnts = self.get_cnts_levels_depths_recs(goobjs)
        self._write_summary_cnts(cnts)

    def _prttex_summary_cnts(self, prt, cnts):
        """Write summary of level and depth counts for active GO Terms."""
        # Count level(shortest path to root) and depth(longest path to root)
        # values for all unique GO Terms.
        prt.write("\n\n% LaTeX Table for GO counts at each level and depth in the GO DAG\n")
        prt.write(r"\begin{table}[bt]" "\n")
        prt.write(r"\begin{tabular}{|r |r |r |r |r |r |r|}" "\n")

        title = self.title.replace('_', r'\_')
        prt.write(r"\hline" "\n")
        prt.write(r"\rowcolor{gray!10}" "\n")
        prt.write(" ".join([r"\multicolumn{7}{|l|}{", title, r"} \\", "\n"]))
        prt.write(r"\hline" "\n")
        prt.write(r"\rowcolor{gray!10}" "\n")
        prt.write(r"Depth &" "\n")
        prt.write(r"\multicolumn{3}{c|}{Depth} &" "\n")
        prt.write(r"\multicolumn{3}{c|}{Level} \\" "\n")
        prt.write(r"\cline{2-7}" "\n")
        prt.write(r"\rowcolor{gray!10}" "\n")
        prt.write(r"or Level & BP & MF & CC & BP & MF & CC \\" "\n")
        prt.write(r"\hline" "\n")

        max_val = max(max(dep for dep in cnts['depth']), max(lev for lev in cnts['level']))
        for i in range(max_val+1):
            vals = ['{:>5}'.format(cnts[desc][i][ns]) for desc in cnts for ns in self.nss]
            self.log.write('{:>02} & {} \\\\\n'.format(i, ' & '.join(vals)))
            if i%2 == 0:
                prt.write(r"\rowcolor{gray!7}" "\n")
        prt.write(r"\hline" "\n")
        prt.write(r"\end{tabular}" "\n")
        prt.write(r"\end{table}" "\n")

    def _write_summary_cnts(self, cnts):
        """Write summary of level and depth counts for active GO Terms."""
        # Count level(shortest path to root) and depth(longest path to root)
        # values for all unique GO Terms.
        max_val = max(max(dep for dep in cnts['depth']),
                      max(lev for lev in cnts['level']))
        self.log.write('Dep <-Depth Counts->  <-Level Counts->\n')
        self.log.write('Lev   BP    MF    CC    BP    MF    CC\n')
        self.log.write('--- ----  ----  ----  ----  ----  ----\n')
        for i in range(max_val+1):
            vals = ['{:>5}'.format(cnts[desc][i][ns]) for desc in sorted(cnts) for ns in self.nss]
            self.log.write('{:>02} {}\n'.format(i, ' '.join(vals)))

    @staticmethod
    def get_cnts_levels_depths_recs(recs):
        """Collect counts of levels and depths in a Group of GO Terms."""
        cnts = cx.defaultdict(lambda: cx.defaultdict(cx.Counter))
        for rec in recs:
            if rec is not None and not rec.is_obsolete:
                cnts['level'][rec.level][rec.namespace] += 1
                cnts['depth'][rec.depth][rec.namespace] += 1
        return cnts

    def get_data(self):
        """Collect counts of GO terms at all levels and depths."""
        # Count level(shortest path to root) and depth(longest path to root)
        # values for all unique GO Terms.
        data = []
        ntobj = cx.namedtuple("NtGoCnt", "Depth_Level BP_D MF_D CC_D BP_L MF_L CC_L")
        cnts = self.get_cnts_levels_depths_recs(set(self.obo.values()))
        max_val = max(max(dep for dep in cnts['depth']), max(lev for lev in cnts['level']))
        for i in range(max_val+1):
            vals = [i] + [cnts[desc][i][ns] for desc in cnts for ns in self.nss]
            data.append(ntobj._make(vals))
        return data

# Copyright (C) 2015-2018, DV Klopfenstein, H Tang, All rights reserved."
