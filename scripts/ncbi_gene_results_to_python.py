#!/usr/bin/env python
"""Print GO terms."""

__copyright__ = "Copyright (C) 2016-2018, DV Klopfenstein, H Tang. All rights reserved."
__author__ = "DV Klopfenstein"

from goatools.cli.ncbi_gene_results_to_python import NCBIgeneToPythonCli


def run():
    """Print GO terms."""
    NCBIgeneToPythonCli().cli()


if __name__ == '__main__':
    run()

# Copyright (C) 2016-2018, DV Klopfenstein, H Tang. All rights reserved.
