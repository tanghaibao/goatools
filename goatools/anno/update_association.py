"""Propagate counts for associations"""

__copyright__ = "Copyright (C) 2016-2019, DV Klopfenstein, H Tang, All rights reserved."
__author__ = "DV Klopfenstein"

import sys
from collections import defaultdict
from goatools.gosubdag.go_tasks import get_go2parents_go2obj
from goatools.anno.broad_gos import NS2GOS_SHORT
from goatools.anno.broad_gos import NS2GOS


def update_association(assc_gene2gos, go2obj, relationships=None, prt=sys.stdout):
    """Add the GO parents of a gene's associated GO IDs to the gene's association."""
    if not assc_gene2gos:
        print('**WARNING: {N} ASSOCATIONS. NO ACTION BY update_association'.format(
            N=len(assc_gene2gos)))
        return
    if prt:
        prt.write("Propagating term counts ")
    # Replaces update_association in GODag
    goids_avail = set(go2obj)
    # Get all assc GO IDs that are current
    assc_goid_sets = assc_gene2gos.values()
    goids_assoc_all = set.union(*assc_goid_sets)
    _chk_goids_notfound(goids_assoc_all, goids_avail)
    # Get the subset of GO objects in the association
    _goids_assoc_cur = goids_assoc_all.intersection(goids_avail)
    _go2obj_assc = {go:go2obj[go] for go in _goids_assoc_cur}
    go2ancestors = get_go2parents_go2obj(_go2obj_assc, relationships, prt)
    # Update the GO sets in assc_gene2gos to include all GO ancestors
    for assc_goids_cur in assc_goid_sets:
        parents = set()
        for goid in assc_goids_cur.intersection(goids_avail):
            if goid in go2ancestors:
                parents.update(go2ancestors[goid])
        assc_goids_cur.update(parents)

def _chk_goids_notfound(goids_assoc_all, goids_avail):
    """Report the number of GO IDs in the association, but not in the GODAG"""
    goids_bad = goids_assoc_all.difference(goids_avail)
    if goids_bad:
        sys.stderr.write("{N} GO IDs NOT FOUND IN ASSOCIATION: {GOs}\n".format(
            N=len(goids_bad), GOs=" ".join(goids_bad)))

def remove_assc_goids(assoc, broad_goids):
    """Remove GO IDs from the association, return a reduced association"""
    actuall_removed_goids = set()
    actuall_removed_genes = set()
    assc_rm = {}
    for geneid, goid_set in assoc.items():
        rm_gos = goid_set.intersection(broad_goids)
        if rm_gos:
            actuall_removed_goids.update(rm_gos)
            reduced_goids = goid_set.difference(rm_gos)
            if reduced_goids:
                assc_rm[geneid] = reduced_goids
            else:
                actuall_removed_genes.add(geneid)
        else:
            assc_rm[geneid] = goid_set
    return {'assoc_reduced':assc_rm,
            'goids_removed':actuall_removed_goids,
            'genes_removed':actuall_removed_genes}


def get_goids_to_remove(goids_or_bool=None):
    """Get GO IDs to remove, with the default being very broad GO IDs"""
    if isinstance(goids_or_bool, (tuple, list, set)):
        return set(goids_or_bool)
    if goids_or_bool is None:
        return set.union(*NS2GOS_SHORT.values())
    if goids_or_bool is True:
        return set.union(*NS2GOS.values())
    return set()


def clean_anno(annots, godag, prt=sys.stdout):
    """Get annotations, gene2gos, for all main GO IDs (not alt) in GO DAG"""
    gene2goset = defaultdict(set)
    go_alts = set()  # For alternate GO IDs
    goids_notfound = set()  # For missing GO IDs
    # Fill go-geneset dict with GO IDs in annotations and their corresponding counts
    for geneid, goids_anno in annots.items():
        # Make a union of all the terms for a gene, if term parents are
        # propagated but they won't get double-counted for the gene
        goids_main = set()
        for goid_anno in goids_anno:
            if goid_anno in godag:
                goid_main = godag[goid_anno].item_id
                if goid_anno != goid_main:
                    go_alts.add(goid_anno)
                goids_main.add(goid_main)
            else:
                goids_notfound.add(goid_anno)
        gene2goset[geneid] = goids_main
    if prt:
        if goids_notfound:
            prt.write("{N} Assc. GO IDs not found in the GODag\n".format(N=len(goids_notfound)))
        prt.write('TermCounts: {N:5} alternate GO IDs in association\n'.format(N=len(go_alts)))
    return dict(gene2goset), go_alts, goids_notfound


# Copyright (C) 2016-2019, DV Klopfenstein, H Tang, All rights reserved.
