"""Utilities for Gene Ontology tasks"""

__copyright__ = "Copyright (C) 2010-2019, DV Klopfenstein et al. All rights reserved."
__author__ = "various"


def get_b2aset(a2bset):
    """Given gene2gos, return go2genes. Given go2genes, return gene2gos."""
    b2aset = {}
    for a_item, bset in a2bset.items():
        for b_item in bset:
            if b_item in b2aset:
                b2aset[b_item].add(a_item)
            else:
                b2aset[b_item] = set([a_item])
    return b2aset


# Copyright (C) 2010-2019, DV Klopfenstein et al. All rights reserved.
