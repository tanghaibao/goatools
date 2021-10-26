"""Utilities for Gene Ontology tasks"""

__copyright__ = "Copyright (C) 2010-present, DV Klopfenstein et al. All rights reserved."
__author__ = "various"

from sys import stdout


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

def read_geneset(fin_genes, prt=stdout):
    """Read a list of genes in a file. Used for study and population genes"""
    genes = set()
    with open(fin_genes) as ifstrm:
        for line in ifstrm:
            line = line.strip()
            if line != '' and line[:1] != '#':
                genes.add(line)
    if genes and next(iter(genes)).isdigit():
        genes = set(int(g) for g in genes)
    if prt:
        prt.write('   {N:6,} READ: {F}\n'.format(N=len(genes), F=fin_genes))
    return genes



# Copyright (C) 2010-present, DV Klopfenstein et al. All rights reserved.
