#!/usr/bin/env python
# -*- coding: UTF-8 -*-

""" Maps gene associations to a 'slim' ontology.

    This roughly implements the functionality of the perl script map2slim.
    See: [http://search.cpan.org/~cmungall/go-perl/scripts/map2slim]

    For a description of GO Slims look here:
    http://geneontology.org/GO.slims.shtml

    For now this does not implement Bucket Terms.
"""

from .obo_parser import GODag


def mapslim(go_term, go_dag, goslim_dag):
    """ Maps a GO term (accession) to it's GO slim terms.

        Parameters:
        - go_term: the accession to be mapped to the slim terms
        - go_dag: the (full) Gene Ontology DAG
        - goslim_dag: the GO Slim DAG

        Returns:
            Two sets:
            direct_ancestors, all_ancestors
        - direct_ancestors: The direct ancestors of the given term that are in
                            the GO Slim. Those are the terms that are not
                            covered by earlier ancestors of the GO Slims in
                            _any_ path (from bottom to top).
        - all_ancestors:    All ancestors of the given term that are part of
                            the GO-Slim terms.

    """
    # check parameters
    if not isinstance(go_dag, GODag):
        raise TypeError("go_dag must be an instance of GODag")
    if not isinstance(goslim_dag, GODag):
        raise TypeError("goslim_dag must be an instance of GODag")
    if go_term not in go_dag:
        raise ValueError("go_term must be an accession that is in the go_dag")

    all_ancestors = set()
    covered_ancestors = set()

    # get all paths for the term in the go_dag
    paths = go_dag.paths_to_top(go_term)
    for path in paths:
        # the next loop needs to run bottom->up, i.e. from the go_term item to
        # the root, thus we need to reverse the list prior to iteration
        path.reverse()

        got_leaf = False
        for term in path:
            if term.id in goslim_dag:
                all_ancestors.add(term.id)
                if got_leaf:
                    covered_ancestors.add(term.id)
                got_leaf = True

    # get the direct ancestors, i.e. those that are not covered by a earlier
    # ancestor of the GO-Slim in _any_ path (in bottom->top order)
    direct_ancestors = all_ancestors - covered_ancestors
    return direct_ancestors, all_ancestors
