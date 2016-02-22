"""Tests downloading and reading of a GO Association File (GAF).

        python test_annotations_gaf.py
"""

__copyright__ = "Copyright (C) 2016, DV Klopfenstein, H Tang. All rights reserved."
__author__ = "DV Klopfenstein"

import os
import sys
from collections import defaultdict
from goatools.associations import read_gaf

def test_gaf_read(log=sys.stdout):
    """Return GO associations from a GAF file. Download if necessary."""
    # Get associations for human(9606), mouse(10090), and fly(7227)
    species_ids = ['goa_human', 'mgi', 'fb']
    # (optional) multi-level dictionary separate associations by taxid
    taxid2asscs = defaultdict(lambda: defaultdict(lambda: defaultdict(set)))
    fin_gafs = dnld_gafs(species_ids)
    for fin_gaf in fin_gafs:
        id2gos = read_gaf(fin_gaf, taxid2asscs=taxid2asscs)
        log.write("  {N} items found in {F}\n".format(N=len(id2gos), F=fin_gaf))
    # Report findings stored in optional taxid dictionary
    for taxid, asscs in taxid2asscs.items():
        num_gene2gos = len(asscs['ID2GOs'])
        num_go2genes = len(asscs['GO2IDs'])
        log.write("{N:>5} GOs and {M:>5} annotated gene ids for tax_id: {TAXID:>6}\n".format(
            TAXID=taxid, N=num_go2genes, M=num_gene2gos))
        # Basic check to ensure gene2go was downloaded and data was returned.
        assert num_gene2gos > 11000
        assert num_go2genes > 6000

def dnld_gafs(species):
    """Download GAF file if necessary."""
    # Example GAF files:
    #   http://geneontology.org/gene-associations/gene_association.goa_human.gz
    #   http://geneontology.org/gene-associations/gene_association.mgi.gz
    #   http://geneontology.org/gene-associations/gene_association.fb.gz
    fin_gafs = []
    for fileid in species:
        gaf_ftp = "http://geneontology.org/gene-associations/gene_association.{S}".format(
            S = fileid)
        gaf_loc = "gene_association.{ID}".format(ID=fileid)
        if not os.path.isfile(gaf_loc):
            os.system("wget {GAF}.gz".format(GAF=gaf))
            os.system("gunzip -f {GAF}.gz".format(GAF=gaf_loc))
        fin_gafs.append(gaf_loc)
    return fin_gafs

if __name__ == '__main__':
    test_gaf_read()

# Copyright (C) 2016, DV Klopfenstein, H Tang. All rights reserved."
