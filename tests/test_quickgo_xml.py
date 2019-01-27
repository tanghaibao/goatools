#!/usr/bin/env python
"""Download a oboxml file from QuickGO and read it."""

from __future__ import print_function

__copyright__ = "Copyright (C) 2015-2018, DV Klopfenstein, H Tang, All rights reserved."
__author__ = "DV Klopfenstein"


from urllib import urlopen
import xmltodict

def test_quickgo():  #### log=sys.stdout):
    """Download a oboxml file from QuickGO and read it."""
    oboxml = urlopen('http://www.ebi.ac.uk/QuickGO/GTerm?id=GO:0003723&format=oboxml')
    print(dir(oboxml))
    print(oboxml.code)
    print(oboxml.headers)
    print(oboxml.info)
    obo_dict = oboxml.read()
    print(obo_dict)
    obo_dict_parse = xmltodict.parse(obo_dict)
    print(obo_dict_parse)


if __name__ == '__main__':
    test_quickgo()

# Copyright (C) 2015-2018, DV Klopfenstein, H Tang, All rights reserved.
