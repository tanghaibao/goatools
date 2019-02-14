"""Semantic Similarity test for Issue #86.

https://github.com/tanghaibao/goatools/issues/86

semantic_similarity & resnik_sim works for few entities but it's giving an error:
    return max(common_parent_go_ids(terms, go), key=lambda t: go[t].depth)
        ValueError: max() arg is an empty sequence

It issues this error when these is no common parent in both provided
entities/genes. Here is one example producing this error
    semantic_similarity(GO:0003676, GO:0007516, godag)

"""

import os
import sys
from goatools.base import get_godag
from goatools.associations import dnld_assc
from goatools.semantic import semantic_distance, semantic_similarity, TermCounts
from goatools.semantic import resnik_sim, lin_sim


def test_top_parent(prt=sys.stdout):
    """Semantic Similarity test for Issue #86."""
    fin_obo = "data/i86.obo"
    branch_dist = 5
    repo = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..")
    godag = get_godag(os.path.join(repo, fin_obo))
    # Get all the annotations from arabidopsis.

    # Calculate the semantic distance and semantic similarity:
    _test_path_same(godag, prt)
    _test_path_parallel(godag, prt)
    _test_path_bp_mf(branch_dist, godag, prt)
    sys.stdout.write("TESTS PASSed: similarity_top_parent\n")

def _test_path_bp_mf(branch_dist, godag, prt):
    """Test distances between BP branch and MF branch."""
    go_mf = 'GO:0003676' # level-03 depth-03 nucleic acid binding [molecular_function]
    go_bp = 'GO:0007516' # level-04 depth-05 hemocyte development [biological_process]
    dst_none = semantic_distance(go_mf, go_bp, godag)
    sim_none = semantic_similarity(go_mf, go_bp, godag)
    assc = dnld_assc("tair.gaf", godag)
    termcounts = TermCounts(godag, assc)
    fmt = '({GO1}, {GO2}) {TYPE:6} score = {VAL}\n'
    sim_r = resnik_sim(go_mf, go_bp, godag, termcounts)
    sim_l = lin_sim(go_mf, go_bp, godag, termcounts)
    if prt is not None:
        prt.write(fmt.format(TYPE='semantic distance', GO1=go_mf, GO2=go_bp, VAL=dst_none))
        prt.write(fmt.format(TYPE='semantic similarity', GO1=go_mf, GO2=go_bp, VAL=sim_none))
        prt.write(fmt.format(TYPE='Resnik similarity', GO1=go_mf, GO2=go_bp, VAL=sim_r))
        prt.write(fmt.format(TYPE='Lin similarity', GO1=go_mf, GO2=go_bp, VAL=sim_l))
    assert dst_none is None
    assert sim_none is None
    assert sim_r is None
    assert sim_l is None
    sim_d = semantic_distance(go_mf, go_bp, godag, branch_dist)
    if prt is not None:
        prt.write(fmt.format(TYPE='semantic distance', GO1=go_mf, GO2=go_bp, VAL=sim_d))
    assert sim_d == godag[go_mf].depth + godag[go_bp].depth + branch_dist

def _test_path_parallel(godag, prt):
    """Test distances between GO IDs on parallel branches."""
    goid_bottom = 'GO:0007516' # BP level-04 depth-05 hemocyte development
    # Test distances up a parallel branch
    goids = [
        'GO:0044763',  # BP level-02 depth-02 single-organism cellular process
        'GO:0008219',  # BP level-03 depth-03 cell death
        'GO:0070997',  # BP level-04 depth-04 neuron death
        'GO:0036475',  # BP level-05 depth-05 neuron death in response to oxidative stress
        'GO:0036476']  # BP level-06 depth-06 neuron death in response to hydrogen peroxide
    fmt = '{DST} semantic_distance between {GO1} and {GO2} on parallel branches\n'
    for dst_exp, goid in enumerate(goids, 3):
        dst_act = semantic_distance(goid_bottom, goid, godag)
        if prt is not None:
            prt.write(fmt.format(DST=dst_act, GO1=goid_bottom, GO2=goid))
        assert dst_act == dst_exp


def _test_path_same(godag, prt):
    """Test distances btween GO IDs on the same path."""
    goid_bottom = 'GO:0007516' # level-04 depth-05 hemocyte development [biological_process]
    # Test distances up the same branch
    goids_bp = [
        'GO:0008150', # level-00 depth-00 biological_process [biological_process]
        'GO:0009987', # level-01 depth-01 cellular process [biological_process]
        'GO:0044763', # level-02 depth-02 single-organism cellular process [biological_process]
        'GO:0048869', # level-03 depth-03 cellular developmental process [biological_process]
        'GO:0048468'] # level-03 depth-04 cell development [biological_process]
    fmt = '{DST} semantic_distance for {GO1} and {GO2} on the same branch\n'
    for dst_exp, goid in enumerate(reversed(goids_bp), 1):
        dst_act = semantic_distance(goid_bottom, goid, godag)
        if prt is not None:
            prt.write(fmt.format(DST=dst_act, GO1=goid_bottom, GO2=goid))
        assert dst_act == dst_exp

if __name__ == '__main__':
    PRT = None if len(sys.argv) != 1 else sys.stdout
    test_top_parent(PRT)
