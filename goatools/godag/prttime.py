"""Print elapsed time"""

__copyright__ = "Copyright (C) 2010-2020, DV Klopfenstein, H Tang, All rights reserved."
__author__ = "DV Klopfenstein"


import os
import sys
import timeit
import datetime
from goatools.obo_parser import GODag
from goatools.base import download_go_basic_obo


# pylint: disable=too-few-public-methods
class GoDagTimed:
    """Load and store GO-DAG. Report elapsed time."""

    def __init__(self, fin_obo, opt_field=None, keep_alt_ids=False):
        self.opt = opt_field  # None causes all fields to read to exp dict
        self.obo = fin_obo
        self._init_dnld_dag()
        self.godag = self.load_dag(self.opt)
        self.go2obj = self.godag if keep_alt_ids else {o.id:o for o in self.godag.values()}

    def _init_dnld_dag(self):
        """If dag does not exist, download it."""
        if not os.path.exists(self.obo):
            download_go_basic_obo(self.obo, loading_bar=None)

    def load_dag(self, opt_fields=None):
        """Run numerous tests for various self.reports."""
        tic = timeit.default_timer()
        dag = GODag(self.obo, opt_fields)
        toc = timeit.default_timer()
        msg = "Elapsed HMS for OBO DAG load: {HMS} OPTIONAL_ATTR({O})\n".format(
            HMS=str(datetime.timedelta(seconds=(toc-tic))), O=opt_fields)
        sys.stdout.write(msg)
        return dag

def prt_hms(tic, msg, prt=sys.stdout):
    """Print elapsed time including Hours, Minutes, and seconds with a user message."""
    hms = str(datetime.timedelta(seconds=(timeit.default_timer()-tic)))
    prt.write("{HMS}: {MSG}\n".format(HMS=hms, MSG=msg))
    return timeit.default_timer()


# Copyright (C) 2010-2020, DV Klopfenstein, H Tang, All rights reserved.
