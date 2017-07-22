"""Semantic Similarity test for Issue #86.

https://github.com/tanghaibao/goatools/issues/86

semantic_similarity & resnik_sim works for few entities but it's giving an error:
    return max(common_parent_go_ids(terms, go), key=lambda t: go[t].depth)
        ValueError: max() arg is an empty sequence

It issues this error when these is no common parent in both provided
entities/genes. Here is one example producing this error
    semantic_similarity(GO:0003676, GO:0007516, godag)

"""

from goatools import obo_parser
from goatools.associations import read_gaf
from goatools.semantic import semantic_similarity

def test_top_parent():
    """Semantic Similarity test for Issue #86."""
    godag = obo_parser.GODag("go-basic.obo")
    # Get all the annotations from arabidopsis.
    associations = read_gaf("http://geneontology.org/gene-associations/gene_association.tair.gz")

    # Calculate the semantic distance and semantic similarity:
    go_id3 = 'GO:0003676'
    go_id4 = 'GO:0007516'
    sim = semantic_similarity(go_id3, go_id4, godag)

if __name__ == '__main__':
    test_top_parent()
