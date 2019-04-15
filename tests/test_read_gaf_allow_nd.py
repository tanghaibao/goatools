#!/usr/bin/env python
"""Read GAF file and allow ND Evidence codes."""

__copyright__ = "Copyright (C) 2016-2019, DV Klopfenstein, H Tang. All rights reserved."
__author__ = "DV Klopfenstein"

import sys
from goatools.associations import read_gaf
from goatools.base import dnld_gaf

def test_gaf_read(log=sys.stdout):
    """Return GO associations from a GAF file. Download if necessary."""
    # On 2017/04/10, there were 3 GO IDs with ND Evidence Codes:
    #
    #    $ cut -f5,7 goa_human.gaf | grep ND | sort | uniq -c
    #        739 GO:0003674      ND
    #        484 GO:0005575      ND
    #        639 GO:0008150      ND

    # Example species_ids: goa_human mgi fb
    fin_gaf = dnld_gaf('goa_human', loading_bar=None)

    # Example 1: Read GAF
    go2ids = read_gaf(fin_gaf, go2geneids=True)
    num_gos_dflt = len(go2ids)
    log.write("Read {N} GOs with all default values\n\n".format(N=num_gos_dflt))

    # Example 2: Read GAF using defaults (No NOT Qualifiers and no ND Evidence Codes)
    go2ids = read_gaf(fin_gaf, go2geneids=True, keep_ND=False, keep_NOT=False)
    log.write("Read {N} GOs; keepif is default in goatools.associations.read_gaf\n\n".format(
        N=len(go2ids)))

    # Example 3: Read GAF allowing GOs with ND Evidence Codes
    go2ids = read_gaf(fin_gaf, go2geneids=True, keep_ND=True)
    log.write("Read {N} GOs; Allow ND Evidence codes\n\n".format(N=len(go2ids)))

    # Example 4: Read GAF allowing all GOs, even those with NOT Qualifiers or ND Evidence Codes
    go2ids = read_gaf(fin_gaf, go2geneids=True, keep_ND=True, keep_NOT=True)
    log.write("Read {N} GOs; Allow ND Evidence codes and NOT Qualifiers\n\n".format(N=len(go2ids)))


if __name__ == '__main__':
    test_gaf_read()

# Copyright (C) 2016-2019, DV Klopfenstein, H Tang. All rights reserved.
