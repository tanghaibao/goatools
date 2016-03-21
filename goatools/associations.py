"""
Routines to read in association file between genes and GO terms.
"""

__copyright__ = "Copyright (C) 2010-2016, H Tang et al. All rights reserved."
__author__ = "various"

from collections import defaultdict
import os
import sys

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

def get_assoc_ncbi_taxids(taxids, force_dnld=False, **kws):
    """Download NCBI's gene2go. Return annotations for user-specified taxid(s)."""
    fin = kws['gene2go'] if 'gene2go' in kws else os.path.join(os.getcwd(), "gene2go")
    dnld_ncbi_gene_file(fin, force_dnld)
    return read_ncbi_gene2go(fin, taxids, **kws)

def dnld_ncbi_gene_file(fin, force_dnld=False, log=sys.stdout):
    """Download a file from NCBI Gene's ftp server."""
    if not os.path.exists(fin) or force_dnld:
        import wget
        import gzip
        fin_dir, fin_base = os.path.split(fin)
        fin_gz = "{F}.gz".format(F=fin_base)
        fin_gz = os.path.join(fin_dir, fin_gz)
        if os.path.exists(fin_gz):
            os.remove(fin_gz)
        fin_ftp = "ftp://ftp.ncbi.nlm.nih.gov/gene/DATA/{F}.gz".format(F=fin_base)
        wget.download(fin_ftp)
        with gzip.open(fin_gz, 'rb') as zstrm:
            log.write("\n  READ:  {F}\n".format(F=fin_gz))
            with open(fin, 'w') as ostrm:
                ostrm.write(zstrm.read())
                log.write("  WROTE: {F}\n".format(F=fin))

def read_ncbi_gene2go(fin_gene2go, taxids=None, **kws):
    """Read NCBI's gene2go. Return gene2go data for user-specified taxids."""
    # Written by DV Klopfenstein
    # kws: taxid2asscs evidence_set
    # Simple associations
    id2gos = defaultdict(set)
    # Optional detailed associations split by taxid and having both ID2GOs & GO2IDs
    # e.g., taxid2asscs = defaultdict(lambda: defaultdict(lambda: defaultdict(set))
    taxid2asscs = kws['taxid2asscs'] if 'taxid2asscs' in kws else None
    evs = kws['evidence_set'] if 'evidence_set' in kws else None
    # By default, return id2gos. User can cause go2geneids to be returned by:
    #   >>> read_ncbi_gene2go(..., go2geneids=True
    b_geneid2gos = not kws['go2geneids'] if 'go2geneids' in kws else True
    if taxids is None: # Default taxid is Human
        taxids = [9606]
    with open(fin_gene2go) as ifstrm:
        for line in ifstrm:
            if line[0] != '#': # Line contains data. Not a comment
                line = line.rstrip() # chomp
                flds = line.split('\t')
                if len(flds) >= 5:
                    taxid_curr, geneid, go_id, evidence, qualifier = line.split('\t')[:5]
                    taxid_curr = int(taxid_curr)
                    # NOT: Used when gene is expected to have function F, but does NOT.
                    # ND : GO function not seen after exhaustive annotation attempts to the gene.
                    if taxid_curr in taxids and qualifier != 'NOT' and evidence != 'ND':
                        # Optionaly specify a subset of GOs based on their evidence.
                        if evs is None or evidence in evs:
                            geneid = int(geneid)
                            if b_geneid2gos:
                                id2gos[geneid].add(go_id)
                            else:
                                id2gos[go_id].add(geneid)
                            if taxid2asscs is not None:
                                taxid2asscs[taxid_curr]['GeneID2GOs'][geneid].add(go_id)
                                taxid2asscs[taxid_curr]['GO2GeneIDs'][go_id].add(geneid)
    return id2gos # return simple associations

def read_gaf(fin_gaf, **kws):
    """Read Gene Association File (GAF). Return data."""
    # Written by DV Klopfenstein
    # kws: taxid2asscs evidence_set
    from goatools.gaf_reader import GafReader
    # Simple associations
    id2gos = defaultdict(set)
    # Optional detailed associations split by taxid and having both ID2GOs & GO2IDs
    taxid2asscs = kws['taxid2asscs'] if 'taxid2asscs' in kws else None
    gafnts = GafReader(fin_gaf).associations
    # Optionaly specify a subset of GOs based on their evidence.
    evs = kws['evidence_set'] if 'evidence_set' in kws else None
    # By default, return id2gos. User can cause go2geneids to be returned by:
    #   >>> read_ncbi_gene2go(..., go2geneids=True
    b_geneid2gos = not kws['go2geneids'] if 'go2geneids' in kws else True
    for ntgaf in gafnts:
        evidence_code = ntgaf.Evidence_Code
        if 'NOT' not in ntgaf.Qualifier and evidence_code != 'ND':
            if evs is None or evidence_code in evs:
                taxid = ntgaf.Taxon[0]
                geneid = ntgaf.DB_ID
                go_id = ntgaf.GO_ID
                if b_geneid2gos:
                    id2gos[geneid].add(go_id)
                else:
                    id2gos[go_id].add(geneid)
                if taxid2asscs is not None:
                    taxid2asscs[taxid]['ID2GOs'][geneid].add(go_id)
                    taxid2asscs[taxid]['GO2IDs'][go_id].add(geneid)
    return id2gos # return simple associations

# Copyright (C) 2010-2016, H Tang et al. All rights reserved."
