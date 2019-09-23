#!/usr/bin/env python
"""Test deprecated location of GoDagTimed"""

import os
import timeit
from goatools.test_data.godag_timed import GoDagTimed
from goatools.test_data.godag_timed import prt_hms
from goatools.base import download_go_basic_obo

REPO = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..")

def test_deprecatedloc_godagtimed():
    """Test deprecated location of GoDagTimed"""
    tic = timeit.default_timer()
    prt_hms(tic, 'prt_hms TESTED')

    fin_go_obo = os.path.join(REPO, "go-basic.obo")
    download_go_basic_obo(fin_go_obo, loading_bar=None)
    GoDagTimed(fin_go_obo)

if __name__ == '__main__':
    test_deprecatedloc_godagtimed()
