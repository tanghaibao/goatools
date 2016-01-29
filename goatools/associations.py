"""
Routines to read in association file between genes and GO terms.
"""

__copyright__ = "Copyright (C) 2016, H Tang. All rights reserved."
__author__ = "various"

from collections import defaultdict
import os

def read_associations(assoc_fn, no_top=False):
    """
    Reads a gene id go term association file. The format of the file
    is as follows:

    AAR1	GO:0005575;GO:0003674;GO:0006970;GO:0006970;GO:0040029
    AAR2	GO:0005575;GO:0003674;GO:0040029;GO:0009845
    ACD5	GO:0005575;GO:0003674;GO:0008219
    ACL1	GO:0005575;GO:0003674;GO:0009965;GO:0010073
    ACL2	GO:0005575;GO:0003674;GO:0009826
    ACL3	GO:0005575;GO:0003674;GO:0009826;GO:0009965

    Also, the following format is accepted (gene ids are repeated):

    AAR1	GO:0005575
    AAR1    GO:0003674
    AAR1    GO:0006970
    AAR2	GO:0005575
    AAR2    GO:0003674
    AAR2    GO:0040029

    :param assoc_fn: file name of the association
    :return: dictionary having keys: gene id, values set of GO terms
    """
    assoc = defaultdict(set)
    top_terms = set(['GO:0008150', 'GO:0003674', 'GO:0005575']) # BP, MF, CC
    for row in open(assoc_fn, 'r'):
        atoms = row.split()
        if len(atoms) == 2:
            gene_id, go_terms = atoms
        elif len(atoms) > 2 and row.count('\t') == 1:
            gene_id, go_terms = row.split("\t")
        else:
            continue
        gos = set(go_terms.split(";"))
        if no_top:
            gos = gos.difference(top_terms)
        assoc[gene_id] |= gos

    return assoc

def get_assoc_ncbi_taxids(taxids, force_dnld=False):
    """Download NCBI's gene2go. Return annotations for user-specified taxid(s)."""
    # Written by DV Klopfenstein, Jan 2016
    import wget
    if not os.path.exists("gene2go") or force_dnld:
        wget.download("ftp://ftp.ncbi.nlm.nih.gov/gene/DATA/gene2go.gz")
        os.system("gunzip gene2go.gz")
    return read_ncbi_gene2go("gene2go", taxids)

def read_ncbi_gene2go(fin_gene2go, taxids):
    """Read NCBI's gene2go. Return gene2go data for user-specified taxids."""
    # Written by DV Klopfenstein, Jan 2016
    taxid2asscs = defaultdict(lambda: defaultdict(lambda: defaultdict(set)))
    with open(fin_gene2go) as ifstrm:
        for line in ifstrm:
            if line[0] != '#': # Line contains data. Not a comment
                line = line.rstrip() # chomp
                flds = line.split('\t')
                if len(flds) >= 5:
                    taxid, geneid, go_id, evidence, qualifier = line.split('\t')[:5]
                    taxid = int(taxid)
                    # NOT: Used when gene is expected to have function F, but does NOT.
                    # ND : GO function not seen after exhaustive annotation attempts to the gene.
                    if taxid in taxids and qualifier != 'NOT' and evidence != 'ND':
                        geneid = int(geneid)
                        taxid2asscs[taxid]['GeneID2GOs'][geneid].add(go_id)
                        taxid2asscs[taxid]['GO2GeneIDs'][go_id].add(geneid)
    return taxid2asscs


