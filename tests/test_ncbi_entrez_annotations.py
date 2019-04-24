#!/usr/bin/env python
"""Tests downloading and reading of the GO annotation file from NCBI Gene.

        python test_NCBI_Entrez_annotations.py
"""

__copyright__ = "Copyright (C) 2016-2019, DV Klopfenstein, H Tang. All rights reserved."
__author__ = "DV Klopfenstein"

import os
import sys
from collections import defaultdict
from goatools.associations import dnld_ncbi_gene_file
from goatools.associations import read_ncbi_gene2go
from goatools.test_data.genes_NCBI_9606_ProteinCoding import GENEID2NT as GeneID2nt_hsa
from goatools.test_data.genes_NCBI_7227_ProteinCoding import GENEID2NT as GeneID2nt_dme

REPO = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..")


def test_ncbi_gene2go(log=sys.stdout):
    """Return GO associations to Entrez GeneIDs. Download if necessary.

       Example report generated with Feb 22, 2013 download of:
         NCBI Gene tables and associations in gene2go

            49672 items found in gene2go from NCBI's ftp server

            taxid    GOs GeneIDs  Description
            ----- ------ -------  -----------
            10090 16,807  18,971  all DNA items
             7227  7,022  12,019  all DNA items
             7227  6,956  10,590  76% GO coverage of 13,919 protein-coding genes
             9606 16,299  18,680  all DNA items
             9606 16,296  18,253  87% GO coverage of 20,913 protein-coding genes

    """
    # Get associations for human(9606), mouse(10090), and fly(7227)
    # (optional) multi-level dictionary separate associations by taxid
    # Simple dictionary containing id2gos
    taxid2asscs = _get_id2gos('gene2go', [9606, 10090, 7227], log)
    taxid2pc = {9606:GeneID2nt_hsa, 7227:GeneID2nt_dme}
    # Report findings
    log.write("   taxid    GOs GeneIDs  Description\n")
    log.write("   ----- ------ -------  -----------\n")
    assert taxid2asscs
    for taxid, asscs in taxid2asscs.items():
        num_gene2gos_all = len(asscs['ID2GOs'])
        num_go2genes_all = len(asscs['GO2IDs'])
        log.write("  {TAXID:>6} {N:>6,} {M:>7,}  all DNA items\n".format(
            TAXID=taxid, N=num_go2genes_all, M=num_gene2gos_all))
        # Basic check to ensure gene2go was downloaded and data was returned.
        assert num_gene2gos_all > 11000
        assert num_go2genes_all > 6000
        if taxid in taxid2pc.keys():
            rpt_coverage(taxid, asscs, taxid2pc[taxid], log)

def _get_id2gos(file_assc, taxids, log):
    """Return associations."""
    taxid2asscs = defaultdict(lambda: defaultdict(lambda: defaultdict(set)))
    fin = os.path.join(REPO, file_assc)
    dnld_ncbi_gene_file(fin, loading_bar=None)
    id2gos = read_ncbi_gene2go(fin, taxids, taxid2asscs=taxid2asscs)
    log.write("  {N} items found in gene2go from NCBI's ftp server\n".format(N=len(id2gos)))
    return taxid2asscs

def rpt_coverage(taxid, asscs, pc2nt, log):
    """Calculate and report GO coverage on protein-coding genes.

       Example report generated with Feb 22, 2013 download of:
         NCBI Gene tables and associations in gene2go

         taxid    GOs GeneIDs  Description
         ----- ------ -------  -----------
          7227  6,956  10,590  76% GO coverage of 13,919 protein-coding genes
          9606 16,296  18,253  87% GO coverage of 20,913 protein-coding genes

    """
    # List of all protein-coding genes have GO terms associated with them
    geneid2gos = asscs['ID2GOs']
    pcgene_w_gos = set(geneid2gos.keys()).intersection(set(pc2nt.keys()))
    num_pcgene_w_gos = len(pcgene_w_gos)
    num_pc_genes = len(pc2nt)
    perc_cov = 100.0*num_pcgene_w_gos/num_pc_genes
    # Get list of GOs associated with protein-coding genes
    gos_pcgenes = set()
    for geneid in pcgene_w_gos:
        gos_pcgenes |= geneid2gos[geneid]
    txt = "  {TAXID:>6} {N:>6,} {M:>7,}  {COV:2.0f}% GO coverage of {TOT:,} protein-coding genes\n"
    log.write(txt.format(
        TAXID=taxid, N=len(gos_pcgenes), M=num_pcgene_w_gos, COV=perc_cov, TOT=num_pc_genes))


if __name__ == '__main__':
    test_ncbi_gene2go()
