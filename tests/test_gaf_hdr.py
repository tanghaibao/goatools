"""Test returning GO Association File (GAF) header version and date."""

__copyright__ = "Copyright (C) 2016-2017, DV Klopfenstein, H Tang. All rights reserved."
__author__ = "DV Klopfenstein"

import os
import sys
from goatools.base import get_gaf_name
from goatools.associations import get_gaf_hdr

REPO = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..")

def main(prt=sys.stdout):
    """Test returning GO Association File (GAF) header version and date."""
    for species in ['goa_human', 'mgi', 'fb']:
        fin_assc = get_gaf_name(species) # goa_human.gaf gene_association.mgi ...
        gafhdr = get_gaf_hdr(fin_assc)
        prt.write("Header from {GAF}:\n{HDR}\n\n".format(GAF=fin_assc, HDR=gafhdr))
        assert "gaf-version" in gafhdr, \
            "UNEXPECTED HEADER FOUND IN {GAF}".format(GAF=fin_assc)

if __name__ == '__main__':
    main()

# Copyright (C) 2016-2017, DV Klopfenstein, H Tang. All rights reserved.
