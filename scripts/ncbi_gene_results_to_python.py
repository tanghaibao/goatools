#!/usr/bin/env python
"""Print GO terms."""

__copyright__ = "Copyright (C) 2016-present, DV Klopfenstein, H Tang. All rights reserved."
__author__ = "DV Klopfenstein"

from goatools.cli.ncbi_gene_results_to_python import main as cli_main


def run():
    """Print GO terms."""
    cli_main()


if __name__ == '__main__':
    run()

# Copyright (C) 2016-present, DV Klopfenstein, H Tang. All rights reserved.
