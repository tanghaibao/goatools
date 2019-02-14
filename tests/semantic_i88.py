#!/usr/bin/env python
"""Test for TermsCount(Go,Association) not working for human genes - goa_human.gaf #88."""
# https://github.com/tanghaibao/goatools/issues/88

from __future__ import print_function

import os
from goatools import obo_parser
from goatools.gosubdag.gosubdag import GoSubDag
from goatools.associations import dnld_assc
from goatools.semantic import semantic_similarity
from goatools.semantic import TermCounts, get_info_content
from goatools.semantic import deepest_common_ancestor
from goatools.semantic import resnik_sim
from goatools.semantic import lin_sim

REPO = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..")


def test_semantic_i88():
    """Computing basic semantic similarities between GO terms."""
    godag = obo_parser.GODag("go-basic.obo")
    goids = set(go for go, o in godag.items() if go == o.id)
    goids = set(godag.keys())
    # Get all the annotations from arabidopsis.
    fin_gaf = os.path.join(REPO, "tair.gaf")
    # dnld_assc includes read_gaf
    associations = dnld_assc(fin_gaf, godag, prt=None)

    # First get the counts and information content for each GO term.
    termcounts = TermCounts(godag, associations)
    gosubdag = GoSubDag(goids, godag, tcntobj=termcounts)

    # Now we can calculate the semantic distance and semantic similarity, as so:
    #       "The semantic similarity between terms GO:0048364 and GO:0044707 is 0.25.
    go_id3 = 'GO:0048364' # BP level-03 depth-04 root development
    go_id4 = 'GO:0044707' # BP level-02 depth-02 single-multicellular organism process
    go_root = deepest_common_ancestor([go_id3, go_id4], godag)
    sim = semantic_similarity(go_id3, go_id4, godag)
    print('\nThe semantic similarity between terms {GO1} and {GO2} is {VAL}.'.format(
        GO1=go_id3, GO2=go_id4, VAL=sim))
    gosubdag.prt_goids([go_root, go_id3, go_id4])

    # Calculate the information content
    go_id = "GO:0048364"
    infocontent = get_info_content(go_id, termcounts)
    print('\nInformation content ({GO}) = {INFO}\n'.format(GO=go_id, INFO=infocontent))

    # Resnik's similarity measure is defined as the information content of the most
    # informative common ancestor. That is, the most specific common parent-term in
    # the GO. Then we can calculate this as follows:
    #       "Resnik similarity score (GO:0048364, GO:0044707) = 4.0540784252
    sim_r = resnik_sim(go_id3, go_id4, godag, termcounts)
    print('Resnik similarity score ({GO1}, {GO2}) = {VAL}'.format(
        GO1=go_id3, GO2=go_id4, VAL=sim_r))

    # Lin similarity score (GO:0048364, GO:0044707) = -0.607721957763
    sim_l = lin_sim(go_id3, go_id4, godag, termcounts)
    print('Lin similarity score ({GO1}, {GO2}) = {VAL}'.format(
        GO1=go_id3, GO2=go_id4, VAL=sim_l))


if __name__ == '__main__':
    test_semantic_i88()
