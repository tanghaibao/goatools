"""Common checks in test data."""

__copyright__ = "Copyright (C) 2010-2019, DV Klopfenstein, H Tang, All rights reserved."
__author__ = "DV Klopfenstein"

from goatools.gosubdag.gosubdag import GoSubDag
from goatools.gosubdag.plot.gosubdag_plot import GoSubDagPlot

class CheckGOs:
    """Check that the dicts of GO IDs are the same"""

    def __init__(self, name, godag):
        self.name = name
        self.godag = godag

    def chk_a2bset(self, exp_a2bset, act_a2bset):
        """Check that the dicts of GO IDs are the same"""
        assert set(exp_a2bset) == set(act_a2bset)
        for goid, exp_goids in exp_a2bset.items():
            act_goids = act_a2bset[goid]
            if act_goids != exp_goids:
                self._err(goid, exp_goids, act_goids)

    def _err(self, goid, exp_goids, act_goids):
        """Report Mismatch Error: Create plot showing relationships"""
        diff_exp = exp_goids.difference(act_goids)
        diff_act = act_goids.difference(exp_goids)
        self._plt(goid, exp_goids, act_goids, diff_exp, diff_act)
        raise RuntimeError('{GO}: EXP[{E}]<->ACT[{A}] EXP({EXP}) ACT({ACT})'.format(
            GO=goid, E=len(exp_goids), A=len(act_goids), EXP=diff_exp, ACT=diff_act))

    # pylint: disable=too-many-arguments
    def _plt(self, goid, exp_goids, act_goids, diff_exp, diff_act):
        """Plot GO IDs, colored by differences in expected and actual"""
        fout_png = '{NAME}_{GO}.png'.format(NAME=self.name, GO=goid.replace(':', ''))
        go_sources = set.union(exp_goids, act_goids, {goid})
        gosubdag = GoSubDag(go_sources, self.godag, relationships=True)
        go2color = {goid: '#c8ffb0'}       # xkcd light light green
        for go_diff in diff_exp:
            go2color[go_diff] = '#cafffb'  # xkcd light light blue
        for go_diff in diff_act:
            go2color[go_diff] = '#ffd1df'  # xkcd light pink
        goploter = GoSubDagPlot(gosubdag, go2color=go2color)
        goploter.plt_dag(fout_png)


# Copyright (C) 2010-2019, DV Klopfenstein, H Tang, All rights reserved.
