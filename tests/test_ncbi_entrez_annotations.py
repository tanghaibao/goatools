"""Tests downloading and reading of the GO annotation file from NCBI Gene.

        python test_NCBI_Entrez_annotations.py
"""

__copyright__ = "Copyright (C) 2016, DV Klopfenstein, H Tang. All rights reserved."
__author__ = "DV Klopfenstein"

import sys
from goatools.associations import get_assoc_ncbi_taxids
from collections import defaultdict
from goatools.test_data.genes_NCBI_9606_ProteinCoding import GeneID2nt as GeneID2nt_hsa
from goatools.test_data.genes_NCBI_7227_ProteinCoding import GeneID2nt as GeneID2nt_dme

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
    taxid2asscs = defaultdict(lambda: defaultdict(lambda: defaultdict(set)))
    # Simple dictionary containing id2gos
    id2gos = get_assoc_ncbi_taxids(taxids=[9606, 10090, 7227], taxid2asscs=taxid2asscs)
    log.write("  {N} items found in gene2go from NCBI's ftp server\n".format(N=len(id2gos)))
    taxid2pc = {9606:GeneID2nt_hsa, 7227:GeneID2nt_dme}
    # Report findings
    log.write("   taxid    GOs GeneIDs  Description\n")
    log.write("   ----- ------ -------  -----------\n")
    for taxid, asscs in taxid2asscs.items():
        num_gene2gos_all = len(asscs['GeneID2GOs'])
        num_go2genes_all = len(asscs['GO2GeneIDs'])
        log.write("  {TAXID:>6} {N:>6,} {M:>7,}  all DNA items\n".format(
            TAXID=taxid, N=num_go2genes_all, M=num_gene2gos_all))
        # Basic check to ensure gene2go was downloaded and data was returned.
        assert num_gene2gos_all > 11000
        assert num_go2genes_all > 6000
        if taxid in taxid2pc.keys():
            rpt_coverage(taxid, asscs, taxid2pc[taxid], log)

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
    geneid2gos = asscs['GeneID2GOs']
    pcgene_w_gos = set(geneid2gos.keys()).intersection(set(pc2nt.keys()))
    num_pcgene_w_gos = len(pcgene_w_gos)
    num_pc_genes = len(pc2nt)
    perc_cov = 100.0*num_pcgene_w_gos/num_pc_genes
    # Get list of GOs associated with protein-coding genes
    gos_pcgenes = set()
    for geneid in pcgene_w_gos:
        gos_pcgenes |= geneid2gos[geneid]
    log.write("  {TAXID:>6} {N:>6,} {M:>7,}  {COV:2.0f}% GO coverage of {TOT:,} protein-coding genes\n".format(
        TAXID=taxid, N=len(gos_pcgenes), M=num_pcgene_w_gos, COV=perc_cov, TOT=num_pc_genes))
    


if __name__ == '__main__':
    test_ncbi_gene2go()
