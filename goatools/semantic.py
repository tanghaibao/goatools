#!/usr/bin/env python
# -*- coding: UTF-8 -*-

"""
Compute semantic similarities between GO terms. Borrowed from book chapter from
Alex Warwick Vesztrocy and Christophe Dessimoz (thanks). For details, please
check out:
notebooks/semantic_similarity.ipynb
"""

import math
from collections import Counter


class TermCounts(object):
    '''
        TermCounts counts the term counts for each
    '''
    def __init__(self, goid, annots):
        '''
            Initialise the counts and
        '''
        # Backup
        self._go = goid

        # Initialise the counters
        self._counts = Counter()
        self._aspect_counts = Counter()

        # Fill the counters...
        self._count_terms(goid, annots)

    def _count_terms(self, goid, annots):
        '''
            Fills in the counts and overall aspect counts.
        '''
        for terms in annots.values(): # key is 'gene'
            # Make a union of all the terms for a gene, if term parents are
            # propagated but they won't get double-counted for the gene
            allterms = set(terms)
            for go_id in terms:
                allterms |= goid[go_id].get_all_parents()
            for parent in allterms:
                self._counts[parent] += 1

        for go_id, child in self._counts.items():
            # Group by namespace
            namespace = goid[go_id].namespace
            self._aspect_counts[namespace] += child

    def get_count(self, go_id):
        '''
            Returns the count of that GO term observed in the annotations.
        '''
        return self._counts[go_id]

    def get_total_count(self, aspect):
        '''
            Gets the total count that's been precomputed.
        '''
        return self._aspect_counts[aspect]

    def get_term_freq(self, go_id):
        '''
            Returns the frequency at which a particular GO term has
            been observed in the annotations.
        '''
        try:
            namespace = self._go[go_id].namespace
            freq = float(self.get_count(go_id)) / float(self.get_total_count(namespace))
            #print self.get_count(go_id), self.get_total_count(namespace), freq
        except ZeroDivisionError:
            freq = 0

        return freq


def get_info_content(go_id, termcounts):
    '''
        Calculates the information content of a GO term.
    '''
    # Get the observed frequency of the GO term
    freq = termcounts.get_term_freq(go_id)

    # Calculate the information content (i.e., -log("freq of GO term")
    return -1.0 * math.log(freq) if freq else 0


def resnik_sim(go_id1, go_id2, goid, termcounts):
    '''
        Computes Resnik's similarity measure.
    '''
    msca = deepest_common_ancestor([go_id1, go_id2], goid)
    return get_info_content(msca, termcounts)


def lin_sim(go_id1, go_id2, goid, termcounts):
    '''
        Computes Lin's similarity measure.
    '''
    sim_r = resnik_sim(go_id1, go_id2, goid, termcounts)

    return (-2*sim_r)/(get_info_content(go_id1, termcounts) + get_info_content(go_id2, termcounts))


def common_parent_go_ids(terms, goid):
    '''
        This function finds the common ancestors in the GO
        tree of the list of terms in the input.
    '''
    # Find candidates from first
    rec = goid[terms[0]]
    candidates = rec.get_all_parents()
    candidates.update({terms[0]})

    # Find intersection with second to nth term
    for term in terms[1:]:
        rec = goid[term]
        parents = rec.get_all_parents()
        parents.update({term})

        # Find the intersection with the candidates, and update.
        candidates.intersection_update(parents)

    return candidates


def deepest_common_ancestor(terms, goid):
    '''
        This function gets the nearest common ancestor
        using the above function.
        Only returns single most specific - assumes unique exists.
    '''
    # Take the element at maximum depth.
    return max(common_parent_go_ids(terms, goid), key=lambda t: goid[t].depth)


def min_branch_length(go_id1, go_id2, goid):
    '''
        Finds the minimum branch length between two terms in the GO DAG.
    '''
    # First get the deepest common ancestor
    dca = deepest_common_ancestor([go_id1, go_id2], goid)

    # Then get the distance from the DCA to each term
    dca_depth = goid[dca].depth
    depth1 = goid[go_id1].depth - dca_depth
    depth2 = goid[go_id2].depth - dca_depth

    # Return the total distance - i.e., to the deepest common ancestor and back.
    return depth1 + depth2


def semantic_distance(go_id1, go_id2, goid):
    '''
        Finds the semantic distance (minimum number of connecting branches)
        between two GO terms.
    '''
    return min_branch_length(go_id1, go_id2, goid)


def semantic_similarity(go_id1, go_id2, goid):
    '''
        Finds the semantic similarity (inverse of the semantic distance)
        between two GO terms.
    '''
    return 1.0 / float(semantic_distance(go_id1, go_id2, goid))
