"""Tests downloading and reading of a GO Association File (GAF).

        python test_annotations_gaf.py
"""

__copyright__ = "Copyright (C) 2016, DV Klopfenstein, H Tang. All rights reserved."
__author__ = "DV Klopfenstein"

import os
import sys
from collections import defaultdict
from goatools.associations import read_gaf
from goatools.base import dnld_gafs

def test_gaf_read(log=sys.stdout):
    """Return GO associations from a GAF file. Download if necessary."""
    # Get associations for human(9606), mouse(10090), and fly(7227)
    species_ids = ['goa_human', 'mgi', 'fb']
    # (optional) multi-level dictionary separate associations by taxid
    taxid2asscs = defaultdict(lambda: defaultdict(lambda: defaultdict(set)))
    local_dir = os.path.dirname(os.path.abspath(__file__))
    for fin_gaf in dnld_gafs(species_ids):
        fin_gaf = os.path.join(local_dir, fin_gaf)
        log.write("\n")
        id2gos = read_gaf(fin_gaf, taxid2asscs=taxid2asscs)
        if "gene_association.mgi" in fin_gaf:
            _chk_key(id2gos, "MGI:")
        log.write("  {N:>6,} IDs found in {F}\n".format(N=len(id2gos), F=fin_gaf))
        go2ids = read_gaf(fin_gaf, go2geneids=True)
        _chk_key(go2ids, "GO:")
        log.write("  {N:>6,} GOs found in {F}\n".format(N=len(go2ids), F=fin_gaf))
    # Report findings stored in optional taxid dictionary
    log.write("\n")
    for taxid, asscs in taxid2asscs.items():
        num_gene2gos = len(asscs.get('ID2GOs'))
        num_go2genes = len(asscs.get('GO2IDs'))
        log.write("{N:>6,} GOs and {M:>6,} annotated gene ids for tax_id: {TAXID:>6}\n".format(
            TAXID=taxid, N=num_go2genes, M=num_gene2gos))
        # Basic check to ensure gene2go was downloaded and data was returned.
        assert num_gene2gos > 11000
        assert num_go2genes > 6000

def _chk_key(a2bs, pattern):
    """Confirm format of dictionary key."""
    for key in a2bs.keys():
        if pattern in key:
            return
        raise RuntimeError("PATTERN({P}) NOT FOUND IN KEY({K})".format(
            P=pattern, K=key))

if __name__ == '__main__':
    test_gaf_read()

# Copyright (C) 2016, DV Klopfenstein, H Tang. All rights reserved.
