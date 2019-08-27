#!/usr/bin/env python
"""Computing basic semantic similarities between GO terms."""

from __future__ import print_function

import os
import itertools
from goatools.base import get_godag
from goatools.associations import dnld_assc
from goatools.semantic import TermCounts
from goatools.semantic import get_info_content
from goatools.semantic import resnik_sim
from goatools.semantic import lin_sim

REPO = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..")


def test_semantic_similarity():
    """Computing basic semantic similarities between GO terms."""
    goids = [
        "GO:0140101",
        "GO:0140097",
        "GO:0140096",
        "GO:0140098",
        "GO:0015318",
        "GO:0140110",
    ]
    # Get all the annotations from arabidopsis.
    associations = [
        ('human', 'goa_human.gaf'),
        ('yeast', 'sgd.gaf'),
    ]


    godag = get_godag(os.path.join(REPO, "go-basic.obo"), loading_bar=None)
    for species, assc_name in associations:  # Limit test numbers for speed
        print()
        # Get all the annotations for the current species
        fin_assc = os.path.join(REPO, assc_name)
        assc_gene2gos = dnld_assc(fin_assc, godag, namespace='MF', prt=None)
        # Calculate the information content of the single term, GO:0048364
        termcounts = TermCounts(godag, assc_gene2gos)

        # Print information values for each GO term
        for goid in sorted(goids):
            infocontent = get_info_content(goid, termcounts)
            term = godag[goid]
            print('{SPECIES} Information content {INFO:8.6f} {NS} {GO} {NAME}'.format(
                SPECIES=species, GO=goid, INFO=infocontent, NS=term.namespace, NAME=term.name))

        # Print semantic similarities between each pair of GO terms
        print("GO #1      GO #2      Resnik Lin")
        print("---------- ---------- ------ -------")
        for go_a, go_b in itertools.combinations(sorted(goids), 2):
            # Resnik's similarity measure is defined as the information content of the most
            # informative common ancestor. That is, the most specific common parent-term in the GO.
            sim_r = resnik_sim(go_a, go_b, godag, termcounts)
            # Lin similarity score (GO:0048364, GO:0044707) = -0.607721957763
            sim_l = lin_sim(go_a, go_b, godag, termcounts)
            print('{GO1} {GO2} {RESNIK:6.4f} {LIN:7.4f}'.format(
                GO1=go_a, GO2=go_b, RESNIK=sim_r, LIN=sim_l))
            assert sim_r >= 0.0, "FATAL RESNIK SCORE: {S}".format(S=sim_r)
            assert sim_l >= 0.0, "FATAL LIN SCORE: {S}".format(S=sim_l)


if __name__ == '__main__':
    test_semantic_similarity()
