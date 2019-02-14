"""Test returning GO Association File (GAF) header version and date."""

__copyright__ = "Copyright (C) 2016-2019, DV Klopfenstein, H Tang. All rights reserved."
__author__ = "DV Klopfenstein"

import os
import sys
from goatools.associations import get_gaf_hdr

REPO = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..")

def main(prt=sys.stdout):
    """Test returning GO Association File (GAF) header version and date."""
    for species in ['goa_human', 'mgi', 'fb']:
        fin_assc = '{SPECIES}.gaf'.format(SPECIES=species) # goa_human.gaf mgi.gaf ...
        gafhdr = get_gaf_hdr(fin_assc)
        prt.write("HEADER FROM {GAF}:\n{HDR}\n\n".format(GAF=fin_assc, HDR=gafhdr))
        assert "gaf-version" in gafhdr, \
            "UNEXPECTED HEADER FOUND IN {GAF}".format(GAF=fin_assc)

if __name__ == '__main__':
    main()

# Copyright (C) 2016-2019, DV Klopfenstein, H Tang. All rights reserved.
