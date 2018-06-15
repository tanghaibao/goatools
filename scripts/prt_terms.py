#!/usr/bin/env python
"""Print GO terms."""

__copyright__ = "Copyright (C) 2016-2018, DV Klopfenstein, H Tang. All rights reserved."
__author__ = "DV Klopfenstein"

from goatools.cli.prt_terms import PrtGOterms


def run():
    """Print GO terms."""
    PrtGOterms().cli()


if __name__ == '__main__':
    run()

# Copyright (C) 2016-2018, DV Klopfenstein, H Tang. All rights reserved.
