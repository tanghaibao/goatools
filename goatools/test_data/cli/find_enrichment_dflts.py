"""Defaults for find_enrichment parseargs to be used in tests."""

__copyright__ = "Copyright (C) 2016-2019, DV Klopfenstein, H Tang, All rights reserved."
__author__ = "DV Klopfenstein"

import os
import collections as cx


# pylint: disable=too-few-public-methods
class ArgsDict(object):
    """Defaults for find_enrichment parseargs to be used in tests."""

    repo = os.path.join(os.path.dirname(os.path.abspath(__file__)), "../../..")

    def __init__(self):
        self.namespace = {
            'annofmt': None,
            'alpha' : 0.05,
            'compare' : False,
            'filenames' : [
                '{REPO}/data/study'.format(REPO=self.repo),
                '{REPO}/data/population'.format(REPO=self.repo),
                '{REPO}/data/association'.format(REPO=self.repo)],
            'goslim' : '{REPO}/goslim_generic.obo'.format(REPO=self.repo),
            'indent' : False,
            'method' : 'bonferroni,sidak,holm,fdr_bh',
            'min_overlap' : 0.7,
            'no_propagate_counts' : False,
            'obo' : '{REPO}/go-basic.obo'.format(REPO=self.repo),
            'outfile' : '{REPO}/goea.txt'.format(REPO=self.repo),
            'outfile_detail' : None,
            'ns': 'BP,MF,CC',
            'pval' : 0.05,
            'pval_field' : 'uncorrected',
            'pvalcalc' : 'fisher',
            'ratio' : None,
            'relationship': False,
            'relationships': None,
            'sections' : None,
            'ev_inc': None,
            'ev_exc': None,
            # BROAD 'remove_goids': None,
        }
        self.ntobj = cx.namedtuple("Namespace", " ".join(self.namespace.keys()))


# Copyright (C) 2016-2019, DV Klopfenstein, H Tang, All rights reserved.
