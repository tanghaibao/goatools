#!/usr/bin/env python
"""Tests GO depth and level hierarchy reporting."""

__copyright__ = "Copyright (C) 2015-2019, DV Klopfenstein, H Tang, All rights reserved."
__author__ = "DV Klopfenstein"

import os
import sys
from collections import defaultdict

from goatools.base import get_godag
from goatools.rpt.rpt_lev_depth import RptLevDepth
from goatools.associations import get_assoc_ncbi_taxids

def test_write_summary_cnts(log=sys.stdout):
    """Print level/depth summaries for various sets of GO terms."""
    fin_obo = os.path.join(os.getcwd(), "go-basic.obo")
    godag = get_godag(fin_obo, loading_bar=None)
    rptobj = RptLevDepth(godag, log)
    # Report level/depth summary for all GOs in a dag
    log.write("\nSummary for all Ontologies:\n")
    rptobj.write_summary_cnts_all()
    # Report level/depth summary for all GOs in human, fly, and mouse
    taxids = [9606, 7227, 10090]
    # (optional) multi-level dictionary separate associations by taxid
    taxid2asscs = defaultdict(lambda: defaultdict(lambda: defaultdict(set)))
    # Get associations for human fly and mouse
    get_assoc_ncbi_taxids(taxids, taxid2asscs=taxid2asscs, loading_bar=None)
    assert taxid2asscs, 'taxid2asscs EMPTY'
    for taxid, assc in taxid2asscs.items():
        log.write("\nSummary for Ontologies for taxid({T}):\n".format(T=taxid))
        go_ids = assc['GO2IDs'].keys()
        rptobj.write_summary_cnts(go_ids)
        log.write("\nSummary for Ontologies for taxid({T}):\n".format(T=taxid))
        go_objs = [godag.get(goid) for goid in go_ids]
        rptobj.write_summary_cnts_goobjs(go_objs)
    # Print GO depth count table for full GO DAG in LaTeX format
    rptobj.prttex_summary_cnts_all(prt=log)

if __name__ == '__main__':
    test_write_summary_cnts()

# Copyright (C) 2015-2019, DV Klopfenstein, H Tang, All rights reserved.
