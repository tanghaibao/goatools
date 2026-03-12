#!/usr/bin/env python
# -*- coding: UTF-8 -*-
"""
Compute semantic similarities between GO terms. Borrowed from book chapter from
Alex Warwick Vesztrocy and Christophe Dessimoz (thanks). For details, please
check out:
notebooks/semantic_similarity.ipynb
"""

import sys
from collections import Counter, defaultdict, deque

from .anno.update_association import clean_anno
from .godag.consts import NAMESPACE2GO, NAMESPACE2NS
from .godag.go_tasks import get_go2ancestors
from .godag.relationship_combos import RelationshipCombos
from .gosubdag.gosubdag import GoSubDag
from .utils import get_b2aset


class TermCounts:
    """
    TermCounts counts the term counts for each
    """

    # pylint: disable=too-many-instance-attributes
    def __init__(self, go2obj, annots, relationships=None, **kws):
        """
        Initialise the counts and
        """
        _prt = kws.get("prt")
        # Handle boolean prt parameter by converting to sys.stdout if True
        if _prt is True:
            _prt = sys.stdout
        elif _prt is False:
            _prt = None
        # Backup
        self.go2obj = go2obj  # Full GODag
        self.annots, go_alts = clean_anno(annots, go2obj, _prt)[:2]
        # Genes annotated to all associated GO, including inherited up ancestors'
        _relationship_set = RelationshipCombos(go2obj).get_set(relationships)
        self.go2genes = self._init_go2genes(_relationship_set, go2obj)
        self.gene2gos = get_b2aset(self.go2genes)
        # Annotation main GO IDs (prefer main id to alt_id)
        self.goids = set(self.go2genes.keys())
        self.gocnts = Counter(
            {go: len(geneset) for go, geneset in self.go2genes.items()}
        )
        # Get total count for each branch: BP MF CC
        self.aspect_counts = {
            "biological_process": self.gocnts.get(
                NAMESPACE2GO["biological_process"], 0
            ),
            "molecular_function": self.gocnts.get(
                NAMESPACE2GO["molecular_function"], 0
            ),
            "cellular_component": self.gocnts.get(
                NAMESPACE2GO["cellular_component"], 0
            ),
        }
        self._init_add_goid_alt(go_alts)
        self.gosubdag = GoSubDag(
            set(self.gocnts.keys()),
            go2obj,
            tcntobj=self,
            relationships=_relationship_set,
            prt=None,
        )
        if _prt:
            self.prt_objdesc(_prt)

    def get_annotations_reversed(self):
        """Return go2geneset for all GO IDs explicitly annotated to a gene"""
        go2genes = get_b2aset(self.annots)
        if go2genes:
            return set.union(*go2genes.values())
        return set()

    def _init_go2genes(self, relationship_set, godag):
        """
        Fills in the genes annotated to each GO, including ancestors

        Due to the ontology structure, gene products annotated to
        a GO Terma are also annotated to all ancestors.
        """
        go2geneset = defaultdict(set)
        go2up = get_go2ancestors(set(godag.values()), relationship_set)
        # Fill go-geneset dict with GO IDs in annotations and their corresponding counts
        for geneid, goids_anno in self.annots.items():
            # Make a union of all the terms for a gene, if term parents are
            # propagated but they won't get double-counted for the gene
            allterms = set()
            for goid_main in goids_anno:
                allterms.add(goid_main)
                if goid_main in go2up:
                    allterms.update(go2up[goid_main])
            # Add 1 for each GO annotated to this gene product
            for ancestor in allterms:
                go2geneset[ancestor].add(geneid)
        return dict(go2geneset)

    def _init_add_goid_alt(self, not_main):
        """
        Add alternate GO IDs to term counts. Report GO IDs not found in GO DAG.
        """
        if not not_main:
            return
        for go_id in not_main:
            if go_id in self.go2obj:
                goid_main = self.go2obj[go_id].item_id
                self.gocnts[go_id] = self.gocnts[goid_main]
                self.go2genes[go_id] = self.go2genes[goid_main]

    def get_count(self, go_id):
        """
        Returns the count of that GO term observed in the annotations.
        """
        return self.gocnts[go_id]

    def get_total_count(self, aspect):
        """
        Gets the total count that's been precomputed.
        """
        return self.aspect_counts[aspect]

    def get_term_freq(self, go_id):
        """
        Returns the frequency at which a particular GO term has
        been observed in the annotations.
        """
        num_ns = float(self.get_total_count(self.go2obj[go_id].namespace))
        return float(self.get_count(go_id)) / num_ns if num_ns != 0 else 0

    def get_gosubdag_all(self, prt=sys.stdout):
        """
        Get GO DAG subset include descendants which are not included in the annotations
        """
        goids = set()
        for gos in self.gosubdag.rcntobj.go2descendants.values():
            goids.update(gos)
        return GoSubDag(
            goids, self.go2obj, self.gosubdag.relationships, tcntobj=self, prt=prt
        )

    def prt_objdesc(self, prt=sys.stdout):
        """Print TermCount object description"""
        ns_tot = sorted(self.aspect_counts.items())
        cnts = [
            "{NS}({N:,})".format(NS=NAMESPACE2NS.get(ns, ns), N=n)
            for ns, n in ns_tot
            if n != 0
        ]
        go_msg = "TermCounts {CNT}".format(CNT=" ".join(cnts))
        prt.write("{GO_MSG} {N:,} genes\n".format(GO_MSG=go_msg, N=len(self.gene2gos)))
        self.gosubdag.prt_objdesc(prt, go_msg)


def get_info_content(go_id, termcounts):
    """
    Retrieve the information content of a GO term.
    """
    if termcounts is None:
        return 0.0
    ntd = termcounts.gosubdag.go2nt.get(go_id)
    return ntd.tinfo if ntd else 0.0


def resnik_sim(go_id1, go_id2, godag, termcounts):
    """
    Computes Resnik's similarity measure.
    """
    goterm1 = godag[go_id1]
    goterm2 = godag[go_id2]
    if goterm1.namespace == goterm2.namespace:
        msca_goid = deepest_common_ancestor([go_id1, go_id2], godag)
        return get_info_content(msca_goid, termcounts)
    return None


def lin_sim(goid1, goid2, godag, termcnts, dfltval=None):
    """
    Computes Lin's similarity measure.
    """
    sim_r = resnik_sim(goid1, goid2, godag, termcnts)
    return lin_sim_calc(goid1, goid2, sim_r, termcnts, dfltval)


def lin_sim_calc(goid1, goid2, sim_r, termcnts, dfltval=None):
    """
    Computes Lin's similarity measure using pre-calculated Resnik's similarities.
    """
    # If goid1 and goid2 are in the same namespace
    if sim_r is not None:
        tinfo1 = get_info_content(goid1, termcnts)
        tinfo2 = get_info_content(goid2, termcnts)
        info = tinfo1 + tinfo2
        # Both GO IDs must be annotated
        if tinfo1 != 0.0 and tinfo2 != 0.0 and info != 0:
            return (2 * sim_r) / info
        # Check if they are identical terms - for identical terms with zero info content, return 1.0
        if termcnts.go2obj[goid1].item_id == termcnts.go2obj[goid2].item_id:
            return 1.0
        # The GOs are separated by the root term, so are not similar
        if sim_r == 0.0:
            return 0.0
    return dfltval


def get_freq_msca(go_id1, go_id2, godag, termcounts):
    """
    Retrieve the frequency of the MSCA of two GO terms.
    """
    try:
        goterm1 = godag[go_id1]
        goterm2 = godag[go_id2]
        if goterm1.namespace == goterm2.namespace:
            msca_goid = deepest_common_ancestor([go_id1, go_id2], godag)
            ntd = termcounts.gosubdag.go2nt.get(msca_goid)
            return ntd.tfreq if ntd else 0
        return 0
    except KeyError:
        return 0


def schlicker_sim(goid1, goid2, godag, termcnts, dfltval=None):
    """
    Computes Schlicker's similarity measure.
    """
    sim_r = resnik_sim(goid1, goid2, godag, termcnts)
    tfreq = get_freq_msca(goid1, goid2, godag, termcnts)
    return schlicker_sim_calc(goid1, goid2, sim_r, tfreq, termcnts, dfltval)


def schlicker_sim_calc(goid1, goid2, sim_r, tfreq, termcnts, dfltval=None):
    """
    Computes Schlicker's similarity measure using pre-calculated Resnik's similarities.
    """
    # If goid1 and goid2 are in the same namespace
    if sim_r is not None:
        tinfo1 = get_info_content(goid1, termcnts)
        tinfo2 = get_info_content(goid2, termcnts)
        info = tinfo1 + tinfo2
        # Both GO IDs must be annotated
        if tinfo1 != 0.0 and tinfo2 != 0.0 and info != 0:
            return (2 * sim_r) / (info) * (1 - tfreq)
        if termcnts.go2obj[goid1].item_id == termcnts.go2obj[goid2].item_id:
            return 1.0 - tfreq
        # The GOs are separated by the root term, so are not similar
        if sim_r == 0.0:
            return 0.0
    return dfltval


def common_parent_go_ids(goids, godag):
    """
    This function finds the common ancestors in the GO
    tree of the list of goids in the input.
    """
    # Find main GO ID candidates from first main or alt GO ID
    rec = godag[goids[0]]
    candidates = rec.get_all_parents()
    candidates.update({rec.item_id})

    # Find intersection with second to nth GO ID
    for goid in goids[1:]:
        rec = godag[goid]
        parents = rec.get_all_parents()
        parents.update({rec.item_id})

        # Find the intersection with the candidates, and update.
        candidates.intersection_update(parents)
    return candidates


def deepest_common_ancestor(goterms, godag):
    """
    This function gets the nearest common ancestor
    using the above function.
    Only returns single most specific - assumes unique exists.
    """
    # Take the element at maximum depth.
    return max(common_parent_go_ids(goterms, godag), key=lambda t: godag[t].depth)


def get_deepest_common_ancestors(goids, godag):
    """
    Returns the set of all deepest (most recent) common ancestors for the
    given GO terms.

    Unlike deepest_common_ancestor(), which returns a single term and may
    behave non-deterministically when multiple common ancestors share the
    same maximum depth, this function returns all of them.

    Parameters
    ----------
    goids : list of str
        List of GO IDs (e.g. ['GO:0048364', 'GO:0044707']).
    godag : GODag
        The GO DAG object (e.g. loaded via GODag('go-basic.obo')).

    Returns
    -------
    set of str
        Set of GO IDs that are the deepest common ancestors of all input terms.
        Returns an empty set if no common ancestor exists.

    Example
    -------
    >>> from goatools.obo_parser import GODag
    >>> from goatools.semantic import get_deepest_common_ancestors
    >>> godag = GODag('go-basic.obo')
    >>> mrca_set = get_deepest_common_ancestors(['GO:0048364', 'GO:0044707'], godag)
    """
    candidates = common_parent_go_ids(goids, godag)
    if not candidates:
        return set()
    max_depth = max(godag[t].depth for t in candidates)
    return {t for t in candidates if godag[t].depth == max_depth}


def min_branch_length(go_id1, go_id2, godag, branch_dist):
    """
    Finds the minimum branch length between two terms in the GO DAG.
    """
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

    if branch_dist is not None:
        return goterm1.depth + goterm2.depth + branch_dist
    return None


def shortest_path(go_id1, go_id2, godag):
    """
    Finds the shortest path length between two GO terms in the same namespace
    by performing a BFS upward through the DAG from each term and finding
    the minimum total distance through any common ancestor.

    Unlike min_branch_length() / semantic_distance(), which compute distance
    using the ``depth`` attribute (longest path from root), this function uses
    a proper BFS traversal and therefore gives correct results in DAGs where
    a term is reachable via multiple paths of different lengths.

    Parameters
    ----------
    go_id1 : str
        First GO ID (e.g. 'GO:0048364').
    go_id2 : str
        Second GO ID (e.g. 'GO:0044707').
    godag : GODag
        The GO DAG object (e.g. loaded via GODag('go-basic.obo')).

    Returns
    -------
    int or None
        The shortest path length (number of edges) between the two GO terms
        through their common ancestors.  Returns ``None`` if the terms are in
        different namespaces (and therefore not connected in the DAG).
        Returns 0 if both GO IDs resolve to the same term.

    Example
    -------
    >>> from goatools.obo_parser import GODag
    >>> from goatools.semantic import shortest_path
    >>> godag = GODag('go-basic.obo')
    >>> dist = shortest_path('GO:0048364', 'GO:0044707', godag)
    """
    goterm1 = godag[go_id1]
    goterm2 = godag[go_id2]
    if goterm1.namespace != goterm2.namespace:
        return None

    def _bfs_up(start_id):
        """BFS upward from start_id through is_a parents; returns {go_id: distance}."""
        distances = {start_id: 0}
        queue = deque([start_id])
        while queue:
            curr_id = queue.popleft()
            curr_dist = distances[curr_id]
            for parent in godag[curr_id].parents:
                pid = parent.item_id
                if pid not in distances:
                    distances[pid] = curr_dist + 1
                    queue.append(pid)
        return distances

    dist1 = _bfs_up(goterm1.item_id)
    dist2 = _bfs_up(goterm2.item_id)

    common = dist1.keys() & dist2.keys()
    if not common:
        return None
    return min(dist1[ancestor] + dist2[ancestor] for ancestor in common)


def semantic_distance(go_id1, go_id2, godag, branch_dist=None):
    """
    Finds the semantic distance (minimum number of connecting branches)
    between two GO terms.
    """
    return min_branch_length(go_id1, go_id2, godag, branch_dist)


def semantic_similarity(go_id1, go_id2, godag, branch_dist=None):
    """
    Finds the semantic similarity (inverse of the semantic distance)
    between two GO terms.
    """
    dist = semantic_distance(go_id1, go_id2, godag, branch_dist)
    if dist is not None:
        return 1.0 / float(dist) if dist != 0 else 1.0
    return None


# 1. Schlicker, Andreas et al.
# "A new measure for functional similarity of gene products based on Gene Ontology"
# BMC Bioinformatics (2006)
#
# 2. Yang, Haixuan et al.
# Improving GO semantic similarity measures by exploring the ontology
#   beneath the terms and modelling uncertainty
# Bioinformatics (2012)
