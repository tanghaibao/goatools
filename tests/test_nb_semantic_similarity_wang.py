#!/usr/bin/env python3
"""Run Jupyter notebook: semantic_similarity_wang"""

from goatools.base import get_godag
from goatools.semsim.termwise.wang import SsWang
from goatools.gosubdag.gosubdag import GoSubDag
from goatools.gosubdag.plot.gosubdag_plot import GoSubDagPlot


# pylint: disable=invalid-name
def test_nb_ss_wang():
    """Run Jupyter notebook: semantic_similarity_wang"""

    go_a = 'GO:0007608'
    go_b = 'GO:0050911'
    go_c = 'GO:0042221'
    relationships = {'part_of'}

    go2txt = {go_a:'GO TERM A', go_b:'GO TERM B', go_c:'GO TERM C'}
    goids = {go_a, go_b, go_c}
    run = Run()

    # pylint: disable=bad-whitespace
    rel_list = [
        (relationships, {}),
        (relationships, {'part_of':.9, 'is_a':.9}),
        ({},            {}),
    ]
    rel_vals = []
    for rels, edge2weight in rel_list:
        wang = SsWang(goids, run.godag, rels, edge2weight)
        wang.prt_cfg()
        fout_png = 'smell_r{N}.png'.format(N=len(rels))
        run.plt(fout_png, goids, rels, go2txt)
        vals = [
            run.get_sim(wang, go_a, go_b),
            run.get_sim(wang, go_a, go_c),
            run.get_sim(wang, go_b, go_c)]
        #rel_vals.append[(rels, vals)]
    print(rel_vals)


class Run:
    """Test details"""

    pattern = ('{GOa} {GOa_name}\n'
               '{GOb} {GOb_name}\n'
               '{VAL:.8f} Wang semantic similarity\n')

    def __init__(self):
        self.godag = get_godag("go-basic.obo", optional_attrs={'relationship'})

    def get_sim(self, wang, go_a, go_b):
        """Calculate Wang's semantic similarity and print details"""
        val = wang.get_sim(go_a, go_b)
        self._print_details(go_a, go_b, val)
        return {'go_a':go_a, 'go_b':go_b, 'val':val}

    def _print_details(self, go_a, go_b, val):
        """Print GO terms and their semantic similarity"""
        print(self.pattern.format(
            GOa=go_a, GOa_name=self.godag[go_a].name,
            GOb=go_b, GOb_name=self.godag[go_b].name,
            VAL=val))

    def plt(self, fout_png, goids, relationships, go2txt):
        """Plot GO terms above resercher GO terms"""
        gosubdag = GoSubDag(goids, self.godag, relationships)
        GoSubDagPlot(gosubdag, go2txt=go2txt).plt_dag(fout_png)


if __name__ == '__main__':
    test_nb_ss_wang()
