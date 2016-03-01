"""Tests GO depth and level hierarchy reporting."""

__copyright__ = "Copyright (C) 2015-2016, DV Klopfenstein, H Tang, All rights reserved."
__author__ = "DV Klopfenstein"

import os
import wget
import sys
from collections import defaultdict

from goatools.obo_parser import GODag
from goatools.rpt_lev_depth import RptLevDepth
from goatools.associations import get_assoc_ncbi_taxids

def test_write_summary_cnts(log=sys.stdout):
    """Print level/depth summaries for various sets of GO terms."""
    obodag = _get_obodag()
    rptobj = RptLevDepth(obodag, log)
    # Report level/depth summary for all GOs in a dag
    log.write("\nSummary for all Ontologies:\n")
    rptobj.write_summary_cnts_all()
    # Report level/depth summary for all GOs in human, fly, and mouse
    taxids = [9606, 7227, 10090]
    # (optional) multi-level dictionary separate associations by taxid
    taxid2asscs = defaultdict(lambda: defaultdict(lambda: defaultdict(set)))
    # Get associations for human fly and mouse
    get_assoc_ncbi_taxids(taxids, taxid2asscs=taxid2asscs)
    for taxid, assc in taxid2asscs.items():
        log.write("\nSummary for Ontologies for taxid({T}):\n".format(T=taxid))
        go_ids = assc['GO2GeneIDs'].keys()
        rptobj.write_summary_cnts(go_ids)
        log.write("\nSummary for Ontologies for taxid({T}):\n".format(T=taxid))
        go_objs = [obodag[goid] for goid in go_ids]
        rptobj.write_summary_cnts_goobjs(go_objs)
    

def _get_obodag():
    """Get GODag from geneontology.org."""
    fin_go_obo = "go-basic.obo"
    if not os.path.exists(fin_go_obo):
        wget.download("http://geneontology.org/ontology/go-basic.obo")
    return GODag(fin_go_obo)

if __name__ == '__main__':
    test_write_summary_cnts()

# Copyright (C) 2015-2016, DV Klopfenstein, H Tang, All rights reserved.
