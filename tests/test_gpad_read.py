#!/usr/bin/env python
"""Test reading GPAD files from Gene Ontology Annotation (GOA) resource http://www.ebi.ac.uk/GOA."""

__copyright__ = "Copyright (C) 2016-2017, DV Klopfenstein, H Tang. All rights reserved."
__author__ = "DV Klopfenstein"

import os
import sys
import collections as cx
from goatools.parsers.dnld_ebi_goa import DnldGoa
from goatools.parsers.gpad_reader import GpadReader

def test_gpad_read(run_full=False, prt=sys.stdout):
    """Test reading GPAD files from GOA source http://www.ebi.ac.uk/GOA."""
    objdnld = DnldGoa()
    species2gpad = _dnld_gpad(objdnld, run_full)
    relations = cx.Counter()
    for org, gpad_file in sorted(species2gpad.items()):
        if org != 'worm': continue
        orgstr = "{ORG} {GPAD}".format(ORG=org, GPAD=os.path.basename(gpad_file))
        prt.write("{GPAD}\n".format(GPAD=orgstr))
        objgpad = GpadReader(gpad_file)
        nts = []
        for ntgpad in objgpad.associations:
            # Assertions are present in the GPAD reader class
            idstr = "{ORG} {ID} {GO}".format(ORG=orgstr, ID=ntgpad.DB_ID, GO=ntgpad.GO_ID)
            if org=="human" and ntgpad.GO_ID == "GO:0004674" and ntgpad.Extension is not None: nts.append(ntgpad)
            if org=="rat" and ntgpad.GO_ID == "GO:0004146" and ntgpad.Extension is not None: nts.append(ntgpad)
            if org=="worm" and ntgpad.GO_ID == "GO:0005634" and ntgpad.Extension is not None: nts.append(ntgpad)
            #prt.write("{E}\n".format(E=ntgpad))
            #prt.write("DB_Reference({E})\n".format(E=ntgpad.DB_Reference))
            #prt.write("ECO_Evidence_Code({E})\n".format(E=ntgpad.ECO_Evidence_Code))
            #if ntgpad.With_From:
            #    prt.write("With_From({E})\n".format(E=ntgpad.With_From))
            #if ntgpad.Taxon:
            #    prt.write("{ID} Taxon({E})\n".format(ID=idstr, E=ntgpad.Taxon))
            #prt.write("Date({E})\n".format(E=ntgpad.Date))
            #prt.write("Assigned_By({E})\n".format(E=ntgpad.Assigned_By))
            if ntgpad.Extension:
                #prt.write("{ID} Extension({E})\n".format(ID=idstr, E=ntgpad.Extension))
                relations += ntgpad.Extension.get_relations_cnt()
            #if ntgpad.Properties:
            #    prt.write("{ID} Properties({E})\n".format(ID=idstr, E=ntgpad.Properties))
        num_ext = len([nt for nt in objgpad.associations if nt.Extension is not None])
        prt.write("{N:8,} of {M:8,} {P:5.2f}% associations have Extensions in {ORG}\n".format(
            N=num_ext, M=objgpad.qty, P=100.*num_ext/objgpad.qty, ORG=org))
        for rel, cnt in objgpad.get_relation_cnt().most_common():
            prt.write("{C:6,} {R}\n".format(C=cnt, R=rel))
    prt.write("\nFULL LIST OF {N} RELATIONS:\n".format(N=len(relations)))
    for rel, cnt in relations.most_common():
        prt.write("{C:10,} {R}\n".format(C=cnt, R=rel))
    # Report GO IDs having extensions of interest.
    for nta in sorted(nts, key=lambda nt: nt.Date, reverse=True):
        txt = "{DB} {GO_ID} {DB_Reference} {Date} {Assigned_By:16}".format(**nta._asdict())
        prt.write("{Ev} {INFO} {EXT}\n".format(
            Ev=nta.Properties['go_evidence'], INFO=txt, EXT=nta.Extension))

def _dnld_gpad(obj, run_full):
    """Return list of downloaded files."""
    species2gpad = {}
    for species in obj.species:
        if not run_full and species == 'uniprot':
            continue
        species2gpad[species] = obj.dnld_goa(species, 'gpa', None)
    return species2gpad

if __name__ == '__main__':
    RUN_FULL = True if len(sys.argv) != 1 else False
    test_gpad_read(RUN_FULL)

# Copyright (C) 2016-2017, DV Klopfenstein, H Tang. All rights reserved."
