#!/usr/bin/env python
# -*- coding: UTF-8 -*-
"""
Compute semantic similarities between GO terms. Borrowed from book chapter from
Alex Warwick Vesztrocy and Christophe Dessimoz (thanks). For details, please
check out:
notebooks/semantic_similarity.ipynb
"""

from __future__ import print_function

import math
from collections import Counter


class TermCounts(object):
    '''
        TermCounts counts the term counts for each
    '''
    def __init__(self, go2obj, annots):
        '''
            Initialise the counts and
        '''
        # Backup
        self.go2obj = go2obj

        # Initialise the counters
        self.gocnts = Counter()
        self.aspect_counts = Counter()

        # Fill the counters...
        self._init_termcounts(annots)


    def _init_termcounts(self, annots):
        '''
            Fill aspect_counts. Find alternate GO IDs that may not be on gocnts.
        '''
        self._init_count_terms(annots)
        self._init_add_goid_alt()


    def _init_count_terms(self, annots):
        '''
            Fills in the counts and overall aspect counts.
        '''
        gonotindag = set()
        gocnts = self.gocnts
        go2obj = self.go2obj
        # Fill gocnts with GO IDs in annotations and their corresponding counts
        for terms in annots.values(): # key is 'gene'
            # Make a union of all the terms for a gene, if term parents are
            # propagated but they won't get double-counted for the gene
            allterms = set()
            for go_id in terms:
                goobj = go2obj.get(go_id, None)
                if goobj is not None:
                    allterms.add(go_id)
                    allterms |= goobj.get_all_parents()
                else:
                    gonotindag.add(go_id)
            for parent in allterms:
                gocnts[parent] += 1
        if gonotindag:
            print("{N} Assc. GO IDs not found in the GODag\n".format(N=len(gonotindag)))


    def _init_add_goid_alt(self):
        '''
            Add alternate GO IDs to term counts.
        '''
        # Fill aspect_counts. Find alternate GO IDs that may not be on gocnts
        goid_alts = set()
        go2cnt_add = {}
        aspect_counts = self.aspect_counts
        gocnts = self.gocnts
        go2obj = self.go2obj
        for go_id, cnt in gocnts.items():
            goobj = go2obj[go_id]
            assert cnt, "NO TERM COUNTS FOR {GO}".format(GO=goobj.id)
            # Was the count set using an alternate GO?
            if go_id != goobj.id:
                go2cnt_add[goobj.id] = cnt
            goid_alts |= goobj.alt_ids
            # Group by namespace
            aspect_counts[goobj.namespace] += cnt
        # If alternate GO used to set count, add main GO ID
        for goid, cnt in go2cnt_add.items():
            gocnts[goid] = cnt
        # Add missing alt GO IDs to gocnts
        for alt_goid in goid_alts.difference(gocnts):
            goobj = go2obj[alt_goid]
            cnt = gocnts[goobj.id]
            assert cnt, "NO TERM COUNTS FOR ALT_ID({GOa}) ID({GO}): {NAME}".format(
                GOa=alt_goid, GO=goobj.id, NAME=goobj.name)
            gocnts[alt_goid] = cnt


    def get_count(self, go_id):
        '''
            Returns the count of that GO term observed in the annotations.
        '''
        return self.gocnts[go_id]

    def get_total_count(self, aspect):
        '''
            Gets the total count that's been precomputed.
        '''
        return self.aspect_counts[aspect]

    def get_term_freq(self, go_id):
        '''
            Returns the frequency at which a particular GO term has
            been observed in the annotations.
        '''
        num_ns = float(self.get_total_count(self.go2obj[go_id].namespace))
        return float(self.get_count(go_id))/num_ns if num_ns != 0 else 0


def get_info_content(go_id, termcounts):
    '''
        Calculates the information content of a GO term.
    '''
    # Get the observed frequency of the GO term
    freq = termcounts.get_term_freq(go_id)

    # Calculate the information content (i.e., -log("freq of GO term")
    return -1.0 * math.log(freq) if freq else 0


def resnik_sim(go_id1, go_id2, godag, termcounts):
    '''
        Computes Resnik's similarity measure.
    '''
    goterm1 = godag[go_id1]
    goterm2 = godag[go_id2]
    if goterm1.namespace == goterm2.namespace:
        msca_goid = deepest_common_ancestor([go_id1, go_id2], godag)
        return get_info_content(msca_goid, termcounts)


def lin_sim(goid1, goid2, godag, termcnts):
    '''
        Computes Lin's similarity measure.
    '''
    sim_r = resnik_sim(goid1, goid2, godag, termcnts)
    return lin_sim_calc(goid1, goid2, sim_r, termcnts)


def lin_sim_calc(goid1, goid2, sim_r, termcnts):
    '''
        Computes Lin's similarity measure using pre-calculated Resnik's similarities.
    '''
    if sim_r is not None:
        info = get_info_content(goid1, termcnts) + get_info_content(goid2, termcnts)
        if info != 0:
            return (-2*sim_r)/info


def common_parent_go_ids(goids, godag):
    '''
        This function finds the common ancestors in the GO
        tree of the list of goids in the input.
    '''
    # Find candidates from first
    rec = godag[goids[0]]
    candidates = rec.get_all_parents()
    candidates.update({goids[0]})

    # Find intersection with second to nth goid
    for goid in goids[1:]:
        rec = godag[goid]
        parents = rec.get_all_parents()
        parents.update({goid})

        # Find the intersection with the candidates, and update.
        candidates.intersection_update(parents)
    return candidates


def deepest_common_ancestor(goterms, godag):
    '''
        This function gets the nearest common ancestor
        using the above function.
        Only returns single most specific - assumes unique exists.
    '''
    # Take the element at maximum depth.
    return max(common_parent_go_ids(goterms, godag), key=lambda t: godag[t].depth)


def min_branch_length(go_id1, go_id2, godag, branch_dist):
    '''
        Finds the minimum branch length between two terms in the GO DAG.
    '''
    # First get the deepest common ancestor
    goterm1 = godag[go_id1]
    goterm2 = godag[go_id2]
    if goterm1.namespace == goterm2.namespace:
        dca = deepest_common_ancestor([go_id1, go_id2], godag)

        # Then get the distance from the DCA to each term
        dca_depth = godag[dca].depth
        depth1 = goterm1.depth - dca_depth
        depth2 = goterm2.depth - dca_depth

        # Return the total distance - i.e., to the deepest common ancestor and back.
        return depth1 + depth2

    elif branch_dist is not None:
        return goterm1.depth + goterm2.depth + branch_dist


def semantic_distance(go_id1, go_id2, godag, branch_dist=None):
    '''
        Finds the semantic distance (minimum number of connecting branches)
        between two GO terms.
    '''
    return min_branch_length(go_id1, go_id2, godag, branch_dist)


def semantic_similarity(go_id1, go_id2, godag, branch_dist=None):
    '''
        Finds the semantic similarity (inverse of the semantic distance)
        between two GO terms.
    '''
    dist = semantic_distance(go_id1, go_id2, godag, branch_dist)
    if dist is not None:
        return 1.0 / float(dist)
