"""Test the loading of the optional GO term fields."""

__copyright__ = "Copyright (C) 2010-2020, DV Klopfenstein, H Tang, All rights reserved."
__author__ = "DV Klopfenstein"

import sys
from goatools.godag.prttime import prt_hms as moved_prt_hms
from goatools.godag.prttime import GoDagTimed as moved_GoDagTimed


# pylint: disable=too-few-public-methods,protected-access
class GoDagTimed:
    """Load and store GO-DAG. Report elapsed time."""

    def __init__(self, fin_obo, opt_field=None, keep_alt_ids=False):
        print('DEPRECATED: GoDagTimed MOVED TO goatools.godag.prttime')
        frm = sys._getframe().f_back.f_code
        print('DEPRECATED GoDagTimed CALLED FROM: {PY} BY {FNC}'.format(
            PY=frm.co_filename, FNC=frm.co_name))
        self.newloc = moved_GoDagTimed(fin_obo, opt_field, keep_alt_ids)
        self.go2obj = self.newloc.go2obj

    def load_dag(self, opt_fields=None):
        """Run numerous tests for various self.reports."""
        return self.newloc.load_dag(opt_fields)

def prt_hms(tic, msg, prt=sys.stdout):
    """Print elapsed time including Hours, Minutes, and seconds with a user message."""
    print('DEPRECATED prt_hms: MOVED TO goatools.godag.prttime')
    frm = sys._getframe().f_back.f_code
    print('DEPRECATED prt_hms CALLED FROM: {PY} BY {FNC}'.format(
        PY=frm.co_filename, FNC=frm.co_name))
    return moved_prt_hms(tic, msg, prt)


# Copyright (C) 2010-2020, DV Klopfenstein, H Tang, All rights reserved.
