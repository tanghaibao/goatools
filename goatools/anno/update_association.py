"""Propagate counts for associations"""

__copyright__ = "Copyright (C) 2016-2019, DV Klopfenstein, H Tang, All rights reserved."
__author__ = "DV Klopfenstein"

import sys
from goatools.gosubdag.go_tasks import get_go2parents_go2obj


def update_association(assc_gene2gos, go2obj, relationships=None, prt=sys.stdout):
    """Add the GO parents of a gene's associated GO IDs to the gene's association."""
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
            parents.update(go2ancestors[goid])
        assc_goids_cur.update(parents)

def _chk_goids_notfound(goids_assoc_all, goids_avail):
    """Report the number of GO IDs in the association, but not in the GODAG"""
    goids_bad = goids_assoc_all.difference(goids_avail)
    if goids_bad:
        sys.stderr.write("{N} GO IDs NOT FOUND IN ASSOCIATION: {GOs}\n".format(
            N=len(goids_bad), GOs=" ".join(goids_bad)))


# Copyright (C) 2016-2019, DV Klopfenstein, H Tang, All rights reserved.
