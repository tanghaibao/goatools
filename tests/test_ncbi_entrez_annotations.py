"""Tests downloading and reading of the GO annotation file from NCBI Gene.

        python test_NCBI_Entrez_annotations.py
"""

__copyright__ = "Copyright (C) 2016, DV Klopfenstein, H Tang. All rights reserved."
__author__ = "DV Klopfenstein"

import sys
from goatools.associations import get_assoc_ncbi_taxids

def test_ncbi_gene2go(log=sys.stdout):
    """Return GO associations to Entrez GeneIDs. Download if necessary."""
    # Get associations for human(9606), mouse(10090), and fly(7227)
    taxid2asscs = get_assoc_ncbi_taxids([9606, 10090, 7227])
    # Report findings
    for taxid, asscs in taxid2asscs.items():
        num_gene2gos = len(asscs['GeneID2GOs'])
        num_go2genes = len(asscs['GO2GeneIDs'])
        log.write("{N:>5} GOs and {M:>5} annotated GeneIDs for tax_id: {TAXID:>6}\n".format(
            TAXID=taxid, N=num_go2genes, M=num_gene2gos))
        # Basic check to ensure gene2go was downloaded and data was returned.
        assert num_gene2gos > 11000
        assert num_go2genes > 6000

if __name__ == '__main__':
    test_ncbi_gene2go()
