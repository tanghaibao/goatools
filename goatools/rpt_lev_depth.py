"""Report the level/depth summaries of all GO terms or a subset of GO terms.

   Level is the length of a shortest path to a GO term.
   Depth is the length of a longest path to a GO term.

    Example:

    >>> obodag = GODag("go-basic.obo")
    >>> reporter = RptLevDepth(obodag)
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

__copyright__ = "Copyright (C) 2015-2016, DV Klopfenstein, H Tang, All rights reserved."
__author__ = "DV Klopfenstein"

import sys
import collections as cx

class RptLevDepth(object):
    """Reports a summary of GO term levels in depths."""

    def __init__(self, obodag, log=sys.stdout):
        self.obo = obodag
        self.log = log

    def write_summary_cnts_all(self):
        """Write summary of level and depth counts for all active GO Terms."""
        cnts = self.get_cnts_levels_depths_recs(set(self.obo.values()))
        self._write_summary_cnts(cnts)

    def write_summary_cnts(self, go_ids):
        """Write summary of level and depth counts for specific GO ids."""
        obo = self.obo
        cnts = self.get_cnts_levels_depths_recs([obo[GO] for GO in go_ids])
        self._write_summary_cnts(cnts)

    def write_summary_cnts_goobjs(self, goobjs):
        """Write summary of level and depth counts for active GO Terms."""
        cnts = self.get_cnts_levels_depths_recs(goobjs)
        self._write_summary_cnts(cnts)

    def _write_summary_cnts(self, cnts):
        """Write summary of level and depth counts for active GO Terms."""
        # Added by DV Klopfenstein
        # Count level(shortest path to root) and depth(longest path to root)
        # values for all unique GO Terms.
        max_val = max(max(dep for dep in cnts['depth']),
                      max(lev for lev in cnts['level']))
        nss = ['biological_process', 'molecular_function', 'cellular_component']
        self.log.write('Dep <-Depth Counts->  <-Level Counts->\n')
        self.log.write('Lev   BP    MF    CC    BP    MF    CC\n')
        self.log.write('--- ----  ----  ----  ----  ----  ----\n')
        for i in range(max_val+1):
            vals = ['{:>5}'.format(cnts[desc][i][ns]) for desc in cnts for ns in nss]
            self.log.write('{:>02} {}\n'.format(i, ' '.join(vals)))

    @staticmethod
    def get_cnts_levels_depths_recs(recs):
        """Collect counts of levels and depths in a Group of GO Terms."""
        cnts = cx.defaultdict(lambda: cx.defaultdict(cx.Counter))
        for rec in recs:
            if not rec.is_obsolete:
                cnts['level'][rec.level][rec.namespace] += 1
                cnts['depth'][rec.depth][rec.namespace] += 1
        return cnts

# Copyright (C) 2015-2016, DV Klopfenstein, H Tang, All rights reserved."
