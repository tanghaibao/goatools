#!/usr/bin/env python
"""Create GO-Dag plots."""

__copyright__ = "Copyright (C) 2016-2018, DV Klopfenstein, H Tang. All rights reserved."
__author__ = "DV Klopfenstein"

import sys

from goatools.cli.main import main as goatools_main


def run():
    """Create GO-Dag plots."""
    goatools_main(["go_plot"] + sys.argv[1:])


if __name__ == '__main__':
    run()

# Copyright (C) 2016-2018, DV Klopfenstein, H Tang. All rights reserved.
