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

from goatools.test_data.optional_attrs import OptionalAttrs


def test_optional_attrs():
    """Test loading optional GO term field, 'synonym'."""
    # Summary for all fields in a GO DAG
    opt_attrs = ['relationship']
    obj = OptionalAttrs("go-basic.obo", opt_attrs)

    # RELATIONSHIP: Stored in a dict with values being sets of GO IDs
    for rel in ['part_of', 'regulates', 'negatively_regulates', 'positively_regulates']:
        obj.chk_relationships_rev(rel, prt=None)
        print("PASSED: relationship: {REL}".format(REL=rel))

    obj.prt_summary()
    obj.chk_get_goterms_upper()
    obj.chk_get_goterms_lower()


def test_no_optional_attrs():
    """Test loading DAG with no optional attributes."""
    obj = OptionalAttrs("go-basic.obo", None)
    obj.chk_no_optattrs()


if __name__ == '__main__':
    test_no_optional_attrs()
    test_optional_attrs()

# Copyright (C) 2010-2018, H Tang et al., All rights reserved.
