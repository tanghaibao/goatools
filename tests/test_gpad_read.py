#!/usr/bin/env python
"""Test reading GPAD files from Gene Ontology Annotation (GOA) resource http://www.ebi.ac.uk/GOA."""

__copyright__ = "Copyright (C) 2016-2019, DV Klopfenstein, H Tang. All rights reserved."
__author__ = "DV Klopfenstein"

import os
import sys
import collections as cx
from goatools.anno.dnld_ebi_goa import DnldGoa
from goatools.anno.gpad_reader import GpadReader
from goatools.base import get_godag


# pylint: disable=too-many-locals
def test_gpad_read(run_desc="mouse", prt=sys.stdout):
    """Test reading GPAD files from GOA source http://www.ebi.ac.uk/GOA."""
    objdnld = DnldGoa()
    species2gpad = _dnld_gpad(objdnld, run_desc)
    # Count Annotation Extension Relations across all species
    relations = cx.Counter()
    godag = get_godag()
    pat = "{N:8,} of {M:8,} {P:5.2f}% associations have Annotation Extensions in {ORG}\n"
    for org, gpad_file in sorted(species2gpad.items()):
        orgstr = "{ORG} {GPAD}".format(ORG=org, GPAD=os.path.basename(gpad_file))
        prt.write("\n{GPAD}\n".format(GPAD=orgstr))
        objgpad = GpadReader(gpad_file, godag=godag)
        for ntgpad in objgpad.associations:
            # Assertions are present in the GPAD reader class
            if ntgpad.Extension:
                relations += ntgpad.Extension.get_relations_cnt()
        num_ext = len([nt for nt in objgpad.associations if nt.Extension is not None])
        # The Extensions field is new in GPAD
        prt.write(pat.format(N=num_ext, M=objgpad.qty, P=100.*num_ext/objgpad.qty, ORG=org))
        for rel, cnt in objgpad.get_relation_cnt().most_common():
            prt.write("    {C:6,} {R}\n".format(C=cnt, R=rel))

    prt.write("\n{N} Annotation Extensions Relations found among all species:\n".format(
        N=len(relations)))
    for rel, cnt in relations.most_common():
        prt.write("{C:10,} {R}\n".format(C=cnt, R=rel))

def _dnld_gpad(objdnld, run_desc):
    """Return list of downloaded files."""
    species2gpad = {}
    species_cur = set(s for s in objdnld.species)
    # Run one species
    if run_desc in species_cur:
        species_cur = set([run_desc])
    # Uniprot is large, so skip it unless specifically asked to include it
    elif run_desc != "inc_uniprot":
        species_cur.remove('uniprot')
    # Download GPAD files for species
    for species in species_cur:
        species2gpad[species] = objdnld.dnld_goa(species, 'gpa', None)
    return species2gpad


if __name__ == '__main__':
    RUN_DESC = "not_uniprot" if len(sys.argv) == 1 else "inc_uniprot"
    test_gpad_read(RUN_DESC)

# Copyright (C) 2016-2019, DV Klopfenstein, H Tang. All rights reserved."
