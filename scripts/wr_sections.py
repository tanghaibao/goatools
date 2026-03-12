#!/usr/bin/env python
"""Create/edit sections files and view GO grouping results."""

__copyright__ = "Copyright (C) 2016-2018, DV Klopfenstein, H Tang. All rights reserved."
__author__ = "DV Klopfenstein"

import sys

from goatools.cli.main import main as goatools_main


def run():
    """Create/edit sections files and view GO grouping results."""
    goatools_main(["wr_sections"] + sys.argv[1:])


if __name__ == '__main__':
    run()

# Copyright (C) 2016-2018, DV Klopfenstein, H Tang. All rights reserved.
