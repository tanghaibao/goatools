#!/usr/bin/env python
"""Code as found in notebooks/semantic_similarity.ipynb."""

from __future__ import print_function

# Computing basic semantic similarities between GO terms

# Adapted from book chapter written by _Alex Warwick Vesztrocy and Christophe Dessimoz_

# How to compute semantic similarity between GO terms.

# First we need to write a function that calculates the minimum number
# of branches connecting two GO terms.

import os
from goatools.base import get_godag
from goatools.associations import dnld_assc
from goatools.semantic import semantic_similarity
from goatools.semantic import TermCounts
from goatools.semantic import get_info_content
from goatools.semantic import deepest_common_ancestor
from goatools.semantic import resnik_sim
from goatools.semantic import lin_sim
from goatools.godag.consts import NS2GO

REPO = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..")

def test_semantic_similarity():
    """Computing basic semantic similarities between GO terms."""
    godag = get_godag(os.path.join(REPO, "go-basic.obo"), loading_bar=None)
    # Get all the annotations from arabidopsis.
    associations = dnld_assc(os.path.join(REPO, 'tair.gaf'), godag)


    # Now we can calculate the semantic distance and semantic similarity, as so:
    #       "The semantic similarity between terms GO:0048364 and GO:0044707 is 0.25.
    go_id3 = 'GO:0048364' # BP level-03 depth-04 root development
    go_id4 = 'GO:0044707' # BP level-02 depth-02 single-multicellular organism process
    sim = semantic_similarity(go_id3, go_id4, godag)
    print('\nThe semantic similarity between terms {GO1} and {GO2} is {VAL}.'.format(
        GO1=go_id3, GO2=go_id4, VAL=sim))
    print(godag[go_id3])
    print(godag[go_id4])

    # Then we can calculate the information content of the single term, <code>GO:0048364</code>.
    #       "Information content (GO:0048364) = 7.75481392334

    # First get the counts of each GO term.
    termcounts = TermCounts(godag, associations)

    # Calculate the information content
    go_id = "GO:0048364"
    infocontent = get_info_content(go_id, termcounts)
    print('\nInformation content ({GO}) = {INFO}\n'.format(GO=go_id, INFO=infocontent))
    assert infocontent, "FATAL INFORMATION CONTENT"

    # Resnik's similarity measure is defined as the information content of the most
    # informative common ancestor. That is, the most specific common parent-term in
    # the GO. Then we can calculate this as follows:
    #       Resnik similarity score (GO:0048364, GO:0044707) = 0.0 because DCA is BP top
    sim_r = resnik_sim(go_id3, go_id4, godag, termcounts)
    dca = deepest_common_ancestor([go_id3, go_id4], godag)
    assert dca == NS2GO['BP']
    assert sim_r == get_info_content(dca, termcounts)
    assert sim_r == 0.0
    print('Resnik similarity score ({GO1}, {GO2}) = {VAL}'.format(
        GO1=go_id3, GO2=go_id4, VAL=sim_r))

    # Lin similarity score (GO:0048364, GO:0044707) = 0.0 because they are similar through BP top
    sim_l = lin_sim(go_id3, go_id4, godag, termcounts)
    print('Lin similarity score ({GO1}, {GO2}) = {VAL}'.format(GO1=go_id3, GO2=go_id4, VAL=sim_l))
    assert sim_l == 0.0, "FATAL LIN SCORE"

    # 
    go_top_cc = NS2GO['CC']
    sim_r = resnik_sim(go_top_cc, go_top_cc, godag, termcounts)
    assert sim_r == 0.0
    sim_l = lin_sim(go_top_cc, go_top_cc, godag, termcounts)
    assert sim_l == 1.0



if __name__ == '__main__':
    test_semantic_similarity()
