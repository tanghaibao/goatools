"""Code as found in notebooks/semantic_similarity.ipynb."""

# Computing basic semantic similarities between GO terms

# Adapted from book chapter written by _Alex Warwick Vesztrocy and Christophe Dessimoz_

# How to compute semantic similarity between GO terms.

# First we need to write a function that calculates the minimum number of branches connecting two GO terms.

from goatools import obo_parser
from goatools.associations import read_gaf
from goatools.semantic import semantic_similarity
from goatools.semantic import TermCounts, get_info_content
from goatools.semantic import resnik_sim
from goatools.semantic import lin_sim

def test_semantic_similarity():
    """Computing basic semantic similarities between GO terms."""
    godag = obo_parser.GODag("go-basic.obo")
    # Get all the annotations from arabidopsis.
    associations = read_gaf("http://geneontology.org/gene-associations/gene_association.tair.gz")


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

    # Resnik's similarity measure is defined as the information content of the most
    # informative common ancestor. That is, the most specific common parent-term in
    # the GO. Then we can calculate this as follows:
    #       "Resnik similarity score (GO:0048364, GO:0044707) = 4.0540784252
    sim_r = resnik_sim(go_id3, go_id4, godag, termcounts)
    print('Resnik similarity score ({GO1}, {GO2}) = {VAL}'.format(GO1=go_id3, GO2=go_id4, VAL=sim_r))

    # Lin similarity score (GO:0048364, GO:0044707) = -0.607721957763
    sim_l = lin_sim(go_id3, go_id4, godag, termcounts)
    print('Lin similarity score ({GO1}, {GO2}) = {VAL}'.format(GO1=go_id3, GO2=go_id4, VAL=sim_l))


if __name__ == '__main__':
    test_semantic_similarity()
