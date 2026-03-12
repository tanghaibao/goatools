#!/usr/bin/env python3
# -*- coding: UTF-8 -*-
"""
python find_enrichment.py study.file population.file gene-association.file

This program returns P-values for functional enrichment in a cluster of study
genes using Fisher's exact test, and corrected for multiple testing (including
Bonferroni, Holm, Sidak, and false discovery rate).

About significance cutoff:
--alpha: test-wise alpha; for each GO term, what significance level to apply
        (most often you don't need to change this other than 0.05 or 0.01)
--pval: experiment-wise alpha; for the entire experiment, what significance
        level to apply after Bonferroni correction
"""

__copyright__ = "Copyright (C) 2010-present, H Tang et al. All rights reserved."
__author__ = "various"

import sys

from goatools.cli.main import main as goatools_main


def main():
    """Run gene enrichment analysis."""
    goatools_main(["find_enrichment"] + sys.argv[1:])


if __name__ == "__main__":
    main()

# Copyright (C) 2010-present, H Tang et al. All rights reserved.
