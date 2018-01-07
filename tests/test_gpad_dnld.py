#!/usr/bin/env python
"""Test downloading files from Gene Ontology Annotation (GOA) resource http://www.ebi.ac.uk/GOA."""

__copyright__ = "Copyright (C) 2016-2018, DV Klopfenstein, H Tang. All rights reserved."
__author__ = "DV Klopfenstein"

import os
import sys
from goatools.anno.dnld_ebi_goa import DnldGoa

def test_ebi_goa_dnld(run_full=False):
    """Test downloading files from GOA source http://www.ebi.ac.uk/GOA."""
    obj = DnldGoa()
    dnld_files = dnld_goa(obj, run_full, ['gpa']) # 'gpi', 'gaf'
    for fout in dnld_files:
        assert os.path.isfile(fout), "FILE({F}) NOT PROPERLY DOWNLOADED FROM {FTP}".format(
            F=fout, FTP=obj.ftp_pub)

def dnld_goa(obj, run_full, exts):
    """Return list of downloaded files."""
    fouts = []
    for species in obj.species:
        if not run_full and species == 'uniprot':
            continue
        for item in [None]: # None, 'complex', 'isoform', 'rna'
            for ext in exts:
                fouts.append(obj.dnld_goa(species, ext, item))
    return fouts

if __name__ == '__main__':
    RUN_FULL = True if len(sys.argv) != 1 else False
    test_ebi_goa_dnld(RUN_FULL)

# Copyright (C) 2016-2018, DV Klopfenstein, H Tang. All rights reserved."
