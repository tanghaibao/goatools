"""Find human genes related to cell cycle."""

import sys
import os
import re
from collections import defaultdict
from goatools.go_search import GoSearch
from goatools.associations import get_assoc_ncbi_taxids
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
    # Because get_assoc_ncbi_taxids returns id2gos, we will opt to 
    # use the (optional) multi-level dictionary separate associations by taxid
    # taxid2asscs contains both GO2GeneIDs and GeneID2GOs.
    taxid2asscs = defaultdict(lambda: defaultdict(lambda: defaultdict(set)))
    get_assoc_ncbi_taxids([taxid], taxid2asscs=taxid2asscs)
    # Initialize GO-search helper object with obo and annotations(go2items)
    srch = GoSearch(fin_go_obo, go2items=taxid2asscs[taxid]['GO2GeneIDs'])
    # Compile search pattern for 'cell cycle'
    cell_cycle = re.compile(r'cell cycle', flags=re.IGNORECASE)
    # Find ALL GOs that have 'cell cycle'. Store results in file.
    fout_allgos = "cell_cycle_gos.log" # Log the search results
    with open(fout_allgos, "w") as prt:
        # Search for 'cell cycle' in GO terms
        gos_cc_all = srch.get_matching_gos(cell_cycle, prt=prt)
        # Researcher carefully reviews GO results and finds GO:0005764(lysosome)
        # in the results when it should not be because the match was found:
        #     cell cycle-independent
        # Researcher removes 'lysosome' from 'cell cycle' results
        # by removing any GOs matching 'cell cycle-independent'
        cell_cycle_ind = re.compile(r'cell cycle.independent', flags=re.IGNORECASE)
        gos_no_cc = srch.get_matching_gos(cell_cycle_ind, gos=gos_cc_all, prt=prt)
        gos = gos_cc_all.difference(gos_no_cc)
        # Add children GOs of cell cycle GOs
        gos_all = srch.add_children_gos(gos)
        log.write('    FOUND {N:>5} GOs:   {F}\n'.format(N=len(gos), F=fout_allgos))
    # Get Entrez GeneIDs for cell cycle GOs
    geneids = srch.get_items(gos_all)
    prt_genes("cell_cycle_genes.txt", geneids, srch, log)
    return geneids

def prt_genes(fout_genes, geneids, srch, log):
    """Print genes."""
    from genes_NCBI_hsa_All import GeneID2nt
    # Print genes related to cell cycle
    fmtstr = "{GeneID:>9} {Symbol:>16} {description}\n"
    nts = [GeneID2nt[geneid] for geneid in sorted(geneids) if geneid in GeneID2nt]
    with open(fout_genes, 'w') as prt:
        prt_txt(prt, nts, fmtstr)
        log.write("    WROTE {N:>5} genes: {FOUT}\n".format(FOUT=fout_genes, N=len(nts)))

if __name__ == '__main__':
    test_cell_cycle()

# Copyright (C) 2010-2016, DV Klopfenstein, H Tang, All rights reserved.
