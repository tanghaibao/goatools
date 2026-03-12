#!/usr/bin/env python
"""Print GO terms."""

__copyright__ = "Copyright (C) 2016-2018, DV Klopfenstein, H Tang. All rights reserved."
__author__ = "DV Klopfenstein"

import sys

from goatools.cli.main import main as goatools_main


def run():
    """Print GO terms."""
    goatools_main(["prt_terms"] + sys.argv[1:])


if __name__ == '__main__':
    run()

# Copyright (C) 2016-2018, DV Klopfenstein, H Tang. All rights reserved.
