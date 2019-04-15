#!/usr/bin/env python
"""Test GoSearch class with both with and without annotations."""

import os
import sys
from collections import defaultdict
from goatools.base import download_go_basic_obo
from goatools.go_search import GoSearch
from goatools.associations import get_assoc_ncbi_taxids

__copyright__ = "Copyright (C) 2010-2019, DV Klopfenstein, H Tang, All rights reserved."
__author__ = "DV Klopfenstein"

REPO = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..")

def test_gosearch(log=sys.stdout):
    """Test GoSearch class with no annotations."""
    taxids = [9606, 10090]
    # Download ontologies and annotations, if necessary
    fin_go_obo = os.path.join(REPO, "go-basic.obo")
    download_go_basic_obo(fin_go_obo, loading_bar=None)
    # Because get_assoc_ncbi_taxids returns id2gos, we will opt to
    # use the (optional) multi-level dictionary separate associations by taxid
    # taxid2asscs contains both GO2IDs and ID2GOs.
    taxid2asscs = defaultdict(lambda: defaultdict(lambda: defaultdict(set)))
    get_assoc_ncbi_taxids(taxids, taxid2asscs=taxid2asscs, loading_bar=None)

    # Initialize GO-search helper object with obo and annotations(go2items)
    for taxid in taxids:
        obj = GoSearch(fin_go_obo, go2items=taxid2asscs[taxid]['GO2IDs'], log=log)
        assert len(obj.obo_dag) > 40000
    GoSearch(fin_go_obo, dict(), log=log)
    assert len(obj.obo_dag) > 40000
    # GoSearch('go.obo', dict(), log=log)


if __name__ == '__main__':
    test_gosearch()

# Copyright (C) 2010-2019, DV Klopfenstein, H Tang, All rights reserved.
