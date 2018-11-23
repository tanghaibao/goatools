#!/usr/bin/env python
"""Compare two or more sets of GO IDs. Best done using sections."""

__copyright__ = "Copyright (C) 2016-2019, DV Klopfenstein, H Tang. All rights reserved."
__author__ = "DV Klopfenstein"

from goatools.cli.compare_gos import CompareGOs


def run():
    """Compare two or more sets of GO IDs. Best done using sections."""
    obj = CompareGOs()
    obj.cli()


if __name__ == '__main__':
    run()

# Copyright (C) 2016-2019, DV Klopfenstein, H Tang. All rights reserved.
