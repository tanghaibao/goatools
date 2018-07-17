#!/usr/bin/env python
"""ArgsDflt minimum overlap."""

from __future__ import print_function

__copyright__ = "Copyright (C) 2010-2018, DV Klopfenstein, H Tang. All rights reserved."

import os
from goatools.base import get_godag
from goatools.cli.find_enrichment import GoeaCliFnc
from goatools.test_data.cli.find_enrichment_dflts import ArgsDict


REPO = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..")

def test_find_enrichment():
    """Recreate run in run.sh."""
    fin_genes = os.path.join(REPO, "data/study")
    pop = set(_.strip() for _ in open(fin_genes) if _.strip())
    stu_orig = pop
    num_pop = len(pop)
    objtest = ArgsDict()
    get_godag(objtest.namespace['obo'], loading_bar=None)
    for min_overlap in [.25, .50, .75]:
        objtest.namespace['min_overlap'] = min_overlap
        args = objtest.ntobj(**objtest.namespace)
        objcli = GoeaCliFnc(args)
        num_stu_in_pop = int(round(min_overlap*num_pop)) + 10
        study = _get_studygenes(stu_orig, num_stu_in_pop)
        overlap = objcli.get_overlap(study, pop)
        print("{N:3} of {M} ({OL}%) in study in pop".format(
            N=num_stu_in_pop, M=num_pop, OL=100.0*overlap))
        objcli.chk_genes(study, pop)
    print("TEST PASSED")

def _get_studygenes(study_orig, num_stu_in_pop):
    """Get a study set having genes not found in the population."""
    study = set()
    for idx, gene in enumerate(study_orig):
        if idx > num_stu_in_pop:
            gene += 'A'
        study.add(gene)
    return study

if __name__ == '__main__':
    test_find_enrichment()

# Copyright (C) 2010-2018, DV Klopfenstein, H Tang. All rights reserved.
