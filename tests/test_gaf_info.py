"""Test returning GO Association File (GAF) header version and date."""

__copyright__ = "Copyright (C) 2016-2017, DV Klopfenstein, H Tang. All rights reserved."
__author__ = "DV Klopfenstein"

import os
import sys
from goatools.base import get_gaf_name
#from goatools.associations import get_gaf_info

REPO = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..")

def main():
    """Test returning GO Association File (GAF) header version and date."""
    for species in ['goa_human', 'mgi', 'fb']:
        fin_assc = get_gaf_name(species) # goa_human.gaf gene_association.mgi ...
        sys.stdout.write("{}\n".format(fin_assc))

if __name__ == '__main__':
    main()

# Copyright (C) 2016-2017, DV Klopfenstein, H Tang. All rights reserved.
