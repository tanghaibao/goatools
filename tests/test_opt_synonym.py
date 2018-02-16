#!/usr/bin/env python
"""Test the loading of the optional GO term field, 'synonym'.

In go-basic.obo fmt(1.2) rel(2018-02-11) 15,329 out of 47,120 GO Terms have synonym types:
  90444 EXACT
  18574 NARROW
  15037 RELATED
   3562 BROAD

"""

from goatools.test_data.optional_attrs import OptionalAttrs


def test_defn():
    """Test loading optional GO term field, 'synonym'."""
    obj = OptionalAttrs("go-basic.obo", None)
    obj.prt_summary()
    # return

    obj = OptionalAttrs("go-basic.obo", "synonym")
    # Check that all GO IDs that should have relationships do have relationships.
    obj.prt_summary()
    # obj.chk_cnt_go()
    # obj.chk_cnt_rel()


if __name__ == '__main__':
    test_defn()
