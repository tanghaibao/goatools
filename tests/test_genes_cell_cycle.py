"""Find human genes related to cell cycle."""

import sys
import os
import re
# sys.path.insert(0, "..") # Use local version of goatools during test
from goatools.go_search import GoSearch
from goatools.associations import get_assoc_ncbi_taxids
from genes_NCBI_hsa_All import GeneID2nt
from goatools.wr_tbl import prt_txt

__copyright__ = "Copyright (C) 2010-2016, DV Klopfenstein, H Tang, All rights reserved."
__author__ = "DV Klopfenstein"

def test_cell_cycle(log=sys.stdout):
    """Test GOEA with local multipletest correction methods."""
    taxid = 9606 # Human annotations
    # Download ontologies and annotations, if necessary
    fin_go_obo = "go-basic.obo"
    if not os.path.exists(fin_go_obo):
        os.system("wget http://geneontology.org/ontology/go-basic.obo")
    taxid2asscs = get_assoc_ncbi_taxids([taxid])
    # Initialize GO-search helper object
    srch = GoSearch(fin_go_obo, go2items=taxid2asscs[taxid]['GO2GeneIDs'])
    # Compile search pattern for 'cell cycle'
    cell_cycle = re.compile(r'cell cycle', flags=re.IGNORECASE)
    # Find ALL GOs that have 'cell cycle'. Store results in file.
    fout_allgos = "cell_cycle_gos.log"
    with open(fout_allgos, "w") as prt:
        gos = srch.get_matching_gos(cell_cycle, prt)
        log.write("    WROTE {N:>5} gos:   {FOUT}\n".format(FOUT=fout_allgos, N=len(gos)))
    # Researcher carefully reviews GO results and finds GO:0005764(lysosome)
    # in the results when it should not be because the match was found:
    #     cell cycle-independent
    # Researcher removes 'lysosome' from 'cell cycle' results
    gos = gos.difference(['GO:0005764'])
    # Add children GOs of cell cycle GOs
    gos_all = srch.add_children_gos(gos)
    # Get Entrez GeneIDs for cell cycle GOs
    geneids = srch.get_items(gos_all)
    # Print genes related to cell cycle
    fmtstr = "{GeneID:>9} {Symbol:>16} {description}\n"
    nts = [GeneID2nt[geneid] for geneid in sorted(geneids) if geneid in GeneID2nt]
    fout_genes = "cell_cycle_genes.txt"
    with open(fout_genes, 'w') as prt:
        prt_txt(prt, nts, fmtstr)
        log.write("    WROTE {N:>5} genes: {FOUT}\n".format(FOUT=fout_genes, N=len(nts)))


if __name__ == '__main__':
    test_cell_cycle()

# Copyright (C) 2010-2016, DV Klopfenstein, H Tang, All rights reserved.
