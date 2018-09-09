#!/usr/bin/env python
"""Test the loading of all optional GO term attributes.

Maximum number of fields on a GO term:

  required:
      1 name
      1 namespace
     10 is_a
     16 alt_id

  optional scalar:
      1 def
      1 comment

  optional set:
     12 subset    goslim_mouse

  optional special:
      3 relationship
    284 synonym   synonym: "peptidase inhibitor complex" EXACT [GOC:bf, GOC:pr]
    827 xref      Wikipedia:Zygotene


In go-basic.obo fmt(1.2) rel(2018-02-11) 15,329 out of 47,120 GO Terms have synonym types:
  90444 EXACT
  18574 NARROW
  15037 RELATED
   3562 BROAD
"""

from __future__ import print_function

import sys
from goatools.test_data.optional_attrs import OptionalAttrs


def test_optional_attrs():
    """Test loading optional GO term field, 'synonym'."""
    args = set(sys.argv[1:])
    prt = sys.stdout if 'prt' in args else None
    exit_if_warning = 'die' in args
    # Summary for all fields in a GO DAG
    opt_attrs = ['def', 'comment', 'subset', 'synonym', 'xref', 'relationship']
    obj = OptionalAttrs("go-basic.obo", opt_attrs)
    obj.prt_summary()

    # SCALAR: Check optional attributes whose information is stored in string
    for attr in ['defn', 'comment']:
        # Check that all GOTerms have a string; ACTUAL matches EXPECTED if present, else ""
        obj.chk_str(attr)
        print("PASSED COUNT TEST: {ATTR}".format(ATTR=attr))

    # SET: Check optional attributes whose information is stored in a set
    # Check that all GO IDs that should have relationships do have relationships.
    # For each GO ID, check that actual count of a set attr equals expected count
    obj.chk_set('subset')
    print("PASSED: subset")

    # RELATIONSHIP: Stored in a dict with values being sets of GO IDs
    obj.chk_relationships()
    print("PASSED: relationship")

    # SYNONYM: Synonyms are stored in a list of namedtuples
    badnts = obj.chk_synonyms(prt)
    _prt_badnts(badnts, exit_if_warning)
    print("PASSED: synonyms")

    # XREF: Stored in a set
    obj.chk_xref(prt)
    print("PASSED: xrefs")

def _prt_badnts(badnts, exit_if_warning):
    """Print bad namedtuples and potentially die."""
    if badnts:
        for goid, badnt in badnts:
            print("**ERROR: BAD NAMEDTUPLE({GO}, {NT})".format(GO=goid, NT=badnt))
        # This option is used to call attention to bad lines in the obo
        if exit_if_warning:
            assert False, "**FATAL: BAD NAMEDTUPLES"

def test_no_optional_attrs():
    """Test loading DAG with no optional attributes."""
    obj = OptionalAttrs("go-basic.obo", None)
    obj.chk_no_optattrs()
    obj = OptionalAttrs("go-basic.obo", [])
    obj.chk_no_optattrs()
    obj = OptionalAttrs("go-basic.obo", set([]))
    obj.chk_no_optattrs()


if __name__ == '__main__':
    test_no_optional_attrs()
    test_optional_attrs()

# Copyright (C) 2010-2018, H Tang et al., All rights reserved.
