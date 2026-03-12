#!/usr/bin/env python
"""Compare two or more sets of GO IDs. Best done using sections."""

__copyright__ = "Copyright (C) 2016-2019, DV Klopfenstein, H Tang. All rights reserved."
__author__ = "DV Klopfenstein"

import sys

from goatools.cli.main import main as goatools_main


def run():
    """Compare two or more sets of GO IDs. Best done using sections."""
    goatools_main(["compare_gos"] + sys.argv[1:])


if __name__ == '__main__':
    run()

# Copyright (C) 2016-2019, DV Klopfenstein, H Tang. All rights reserved.
