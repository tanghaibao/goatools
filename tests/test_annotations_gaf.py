#!/usr/bin/env python
"""Tests downloading and reading of a GO Association File (GAF).

        python test_annotations_gaf.py
"""

__copyright__ = "Copyright (C) 2016-2019, DV Klopfenstein, H Tang. All rights reserved."
__author__ = "DV Klopfenstein"

import os
import sys

from collections import defaultdict
from goatools.associations import read_gaf
from goatools.base import dnld_gafs

def test_gaf_read(log=sys.stdout):
    """Return GO associations from a GAF file. Download if necessary."""
    # Get associations for human(9606), mouse(10090), and fly(7227)
    # Read GAF associations
    msg = "Read GAF associations; keepif == None (default behavior)"
    species_ids = ['goa_human', 'mgi', 'fb']
    _test_gaf_read(msg, species_ids, None, log)
    # Read GAF associations
    msg = "Read GAF associations; keepif is default in goatools.associations.read_gaf"
    keepif = lambda nt: 'NOT' not in nt.Qualifier and nt.Evidence_Code != 'ND'
    _test_gaf_read(msg, species_ids, keepif, log)
    # Read GAF associations, allowing ND Evidence codes
    msg = "Read GAF associations; Allow ND Evidence codes"
    keepif = lambda nt: 'NOT' not in nt.Qualifier
    _test_gaf_read(msg, species_ids, keepif, log)
    # Read GAF associations, allowing ND entries and NOT Qualifiers
    msg = "Read GAF associations; Allow ND Evidence codes and NOT Qualifiers"
    keepif = lambda nt: True
    #_test_gaf_read(msg, species_ids, keepif, log)
    # Limit number of tests for speed
    _test_gaf_read(msg, species_ids[-1:], keepif, log)

def _test_gaf_read(msg, species_ids, keepif, log=sys.stdout):
    # (optional) multi-level dictionary separate associations by taxid
    taxid2asscs = defaultdict(lambda: defaultdict(lambda: defaultdict(set)))
    local_dir = os.path.dirname(os.path.abspath(__file__))
    for fin_gaf in dnld_gafs(species_ids, loading_bar=None):
        fin_gaf = os.path.join(local_dir, fin_gaf)
        log.write("\n")
        id2gos_bp = read_gaf(fin_gaf, taxid2asscs=taxid2asscs, keepif=keepif)
        id2gos_all = read_gaf(fin_gaf, taxid2asscs=taxid2asscs, keepif=keepif, namespace='all')
        assert len(id2gos_all) > len(id2gos_bp)
        if "mgi.gaf" in fin_gaf:
            _chk_key(id2gos_bp, "MGI:")
        log.write("  {N:>6,} IDs found in BP  {F}\n".format(N=len(id2gos_bp), F=fin_gaf))
        log.write("  {N:>6,} IDs found in ALL {F}\n".format(N=len(id2gos_all), F=fin_gaf))
        go2ids = read_gaf(fin_gaf, go2geneids=True, keepif=keepif)
        _chk_key(go2ids, "GO:")
        log.write("  {N:>6,} GOs found in {F}\n".format(N=len(go2ids), F=fin_gaf))
    # Report findings stored in optional taxid dictionary
    log.write("\n{MSG}\n".format(MSG=msg))
    txtpat = "  {N:>6,} GOs and {M:>6,} annotated gene ids for tax_id: {TAXID:>6}\n"
    for taxid, asscs in taxid2asscs.items():
        num_gene2gos = len(asscs.get('ID2GOs'))
        num_go2genes = len(asscs.get('GO2IDs'))
        log.write(txtpat.format(TAXID=taxid, N=num_go2genes, M=num_gene2gos))
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

# Copyright (C) 2016-2019, DV Klopfenstein, H Tang. All rights reserved.
