"""Find human genes related to cell cycle."""

import sys
import os
import re
from collections import defaultdict
from goatools.base import download_go_basic_obo
from goatools.go_search import GoSearch
from goatools.associations import get_assoc_ncbi_taxids
from goatools.wr_tbl import prt_txt

__copyright__ = "Copyright (C) 2010-2019, DV Klopfenstein, H Tang, All rights reserved."
__author__ = "DV Klopfenstein"

def test_cell_cycle(taxid=9606, log=sys.stdout):
    """Get all genes related to cell cycle. Write results to file."""
    geneids = get_genes_cell_cycle(taxid, log)
    fout = "cell_cycle_genes_{TAXID}.txt".format(TAXID=taxid)
    prt_genes(fout, geneids, taxid, log)

def get_genes_cell_cycle(taxid=9606, log=sys.stdout):
    """Test GOEA with local multipletest correction methods for cell cycle."""
    # Download ontologies and annotations, if necessary
    fin_go_obo = os.path.join(os.getcwd(), "go-basic.obo")
    download_go_basic_obo(fin_go_obo, loading_bar=None)
    # Because get_assoc_ncbi_taxids returns id2gos, we will opt to
    # use the (optional) multi-level dictionary separate associations by taxid
    # taxid2asscs contains both GO2IDs and ID2GOs.
    taxid2asscs = defaultdict(lambda: defaultdict(lambda: defaultdict(set)))
    get_assoc_ncbi_taxids([taxid], taxid2asscs=taxid2asscs, loading_bar=None)

    # Initialize GO-search helper object with obo and annotations(go2items)
    srch = GoSearch(fin_go_obo, go2items=taxid2asscs[taxid]['GO2IDs'])
    # Compile search pattern for 'cell cycle'
    cell_cycle = re.compile(r'cell cycle', flags=re.IGNORECASE)
    # Find ALL GOs that have 'cell cycle'. Store results in file.
    fout_allgos = "cell_cycle_gos_{TAXID}.log".format(TAXID=taxid)
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
        if log is not None:
            log.write('    taxid {TAXID:>5}\n'.format(TAXID=taxid))
            log.write('    FOUND {N:>5} GOs:   {F}\n'.format(
                N=len(gos_all), F=fout_allgos))
    # Get Entrez GeneIDs for cell cycle GOs
    geneids = srch.get_items(gos_all)
    return geneids

def prt_genes(fout_genes, geneids, taxid, log):
    """Print 'cell cycle' geneids, with or without Symbol and description information."""
    fin_symbols = "genes_NCBI_{TAXID}_All.py".format(TAXID=taxid)
    # If gene Symbol information is available, print geneid and Symbol
    if os.path.isfile(fin_symbols):
        import importlib
        module_name = "".join(["goatools.test_data.", fin_symbols[:-3]])
        module = importlib.import_module(module_name)
        geneid2nt = module.GENEID2NT
        fmtstr = "{GeneID:>9} {Symbol:<16} {description}\n"
        nts = [geneid2nt[geneid] for geneid in sorted(geneids) if geneid in geneid2nt]
        with open(fout_genes, 'w') as prt:
            prt_txt(prt, nts, fmtstr)
            if log is not None:
                log.write("    WROTE {N:>5} genes: {FOUT}\n".format(FOUT=fout_genes, N=len(nts)))
    # Just print geneids
    else:
        with open(fout_genes, 'w') as prt:
            for geneid in geneids:
                prt.write("{geneid}\n".format(geneid=geneid))
            if log is not None:
                log.write("    WROTE {N:>5} genes: {FOUT}\n".format(
                    FOUT=fout_genes, N=len(geneids)))

if __name__ == '__main__':
    test_cell_cycle(9606)
    test_cell_cycle(10090)

# Copyright (C) 2010-2019, DV Klopfenstein, H Tang, All rights reserved.
