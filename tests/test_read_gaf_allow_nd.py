"""Read GAF file and allow ND Evidence codes."""

__copyright__ = "Copyright (C) 2016-2017, DV Klopfenstein, H Tang. All rights reserved."
__author__ = "DV Klopfenstein"

import sys
from goatools.associations import read_gaf
from goatools.base import dnld_gaf

def test_gaf_read(log=sys.stdout):
    """Return GO associations from a GAF file. Download if necessary."""
    # Example species_ids: goa_human mgi fb
    fin_gaf = dnld_gaf('goa_human')

    # Example 1: Read GAF
    go2ids = read_gaf(fin_gaf, go2geneids=True)
    log.write("Read {N} GOs with all default values\n\n".format(N=len(go2ids)))

    # Example 2: Read GAF using defaults (No NOT Qualifiers and no ND Evidence Codes)
    keepif = lambda nt: 'NOT' not in nt.Qualifier and nt.Evidence_Code != 'ND'
    go2ids = read_gaf(fin_gaf, go2geneids=True, keepif=keepif)
    log.write("Read {N} GOs; keepif is default in goatools.associations.read_gaf\n\n".format(
        N=len(go2ids)))

    # Example 3: Read GAF allowing GOs with ND Evidence Codes
    keepif = lambda nt: 'NOT' not in nt.Qualifier
    go2ids = read_gaf(fin_gaf, go2geneids=True, keepif=keepif)
    log.write("Read {N} GOs; Allow ND Evidence codes\n\n".format(N=len(go2ids)))

    # Example 4: Read GAF allowing all GOs, even those with NOT Qualifiers or ND Evidence Codes
    keepif = lambda nt: True
    go2ids = read_gaf(fin_gaf, go2geneids=True, keepif=keepif)
    log.write("Read {N} GOs; Allow ND Evidence codes and NOT Qualifiers\n\n".format(N=len(go2ids)))


if __name__ == '__main__':
    test_gaf_read()

# Copyright (C) 2016-2017, DV Klopfenstein, H Tang. All rights reserved.
