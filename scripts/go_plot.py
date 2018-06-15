#!/usr/bin/env python
"""Create GO-Dag plots."""

__copyright__ = "Copyright (C) 2016-2018, DV Klopfenstein, H Tang. All rights reserved."
__author__ = "DV Klopfenstein"

from goatools.cli.gosubdag_plot import PlotCli


def run():
    """Create GO-Dag plots."""
    PlotCli().cli()


if __name__ == '__main__':
    run()

# Copyright (C) 2016-2018, DV Klopfenstein, H Tang. All rights reserved.
