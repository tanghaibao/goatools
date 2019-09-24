#!/usr/bin/env python
"""Test loading of relationships, like part_of, into TermCounts"""

import os
import sys
from goatools.obo_parser import GODag
from goatools.anno.idtogos_reader import IdToGosReader
from goatools.semantic import TermCounts
from goatools.gosubdag.plot.gosubdag_plot import GoSubDagPlot


REPO = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..")

def test_tcntobj_relationships(do_plt=False):
    """Test loading of relationships, like part_of, into TermCounts"""
    # Filenames
    fin_obo = os.path.join(REPO, "tests/data/yangRWC/fig2a.obo")
    fin_anno = os.path.join(REPO, "tests/data/yangRWC/fig2a.anno")
    fout_png_r0 = os.path.join(REPO, 'yang_fig2a_r0.png')
    fout_png_r1 = os.path.join(REPO, 'yang_fig2a_r1.png')
    relationships = {'part_of',}

    # Load ontologies
    go2obj = GODag(fin_obo, optional_attrs=['relationship'])

    # Load annotations
    assoc = IdToGosReader(fin_anno, godag=go2obj).get_id2gos('CC')

    # Count genes annotated to GO terms w and wo/relationships
    tcntobj_r0 = TermCounts(go2obj, assoc)
    # relationship: G (GO:0000007) is part_of F (GO:0000006)
    tcntobj_r1 = TermCounts(go2obj, assoc, relationships)

    # Check results
    # Adding relationships does not change the total count of genes:
    assert tcntobj_r0.gocnts['GO:0005575'] == tcntobj_r1.gocnts['GO:0005575']
    # Counts without relationships:
    assert tcntobj_r0.gocnts['GO:0000002'] == 40  # GO Term B
    assert tcntobj_r0.gocnts['GO:0000006'] == 10  # GO Term F
    # Counts with relationships: F counts G's 30 genes, so does B
    assert tcntobj_r1.gocnts['GO:0000002'] == 70  # GO Term B
    assert tcntobj_r1.gocnts['GO:0000006'] == 40  # GO Term F

    # Optionally visualize the difference between term counts w and wo/relationships
    if do_plt:
        go2txt_r0 = {nt.GO:'tcnt={}'.format(nt.tcnt) for nt in tcntobj_r0.gosubdag.go2nt.values()}
        GoSubDagPlot(tcntobj_r0.gosubdag, go2txt=go2txt_r0).plt_dag(fout_png_r0)
        go2txt_r1 = {nt.GO:'tcnt={}'.format(nt.tcnt) for nt in tcntobj_r1.gosubdag.go2nt.values()}
        GoSubDagPlot(tcntobj_r1.gosubdag, go2txt=go2txt_r1).plt_dag(fout_png_r1)


if __name__ == '__main__':
    test_tcntobj_relationships(len(sys.argv) != 1)
