#!/usr/bin/env python3
"""Test running an enrichment using any annotation file format."""

from __future__ import print_function

__copyright__ = "Copyright (C) 2010-2019, DV Klopfenstein, H Tang. All rights reserved."

import os
import itertools
from goatools.base import get_godag
from goatools.associations import dnld_annofile
from goatools.anno.factory import get_objanno
from goatools.gosubdag.gosubdag import GoSubDag
from goatools.goea.go_enrichment_ns import GOEnrichmentStudyNS

REPO = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..")


def test_find_enrichment():
    """RUn an enrichments using all annotation file formats"""

    godag = get_godag("go-basic.obo", optional_attrs=['relationship'])
    gos = _get_enriched_goids('GO:0006959', godag)  # GO IDs related to humoral response

    # pylint: disable=superfluous-parens
    print('- DOWNLOAD AND LOAD -----------------------------------------------')
    annoobjs = [
        _get_objanno('gene2go', taxid=10090),
        _get_objanno('gene2go', taxid=9606),
        _get_objanno('goa_human.gaf'),
        _get_objanno('goa_human.gpad', godag=godag),
        _get_objanno('data/association', anno_type='id2gos', godag=godag),
    ]

    for obj in annoobjs:
        ns2assc = obj.get_ns2assc()
        pop = list(itertools.chain.from_iterable(ns2assc.values()))
        print('{N:6,} population IDs'.format(N=len(pop)))
        enriched = set(nt.DB_ID for nt in obj.associations if nt.GO_ID in gos)
        objgoeans = _get_objgoeans(pop, ns2assc, godag)
        results = objgoeans.run_study(enriched)
        print('{N} results'.format(N=len(results)))
        # Run one branch
        bp2assc = {'BP': ns2assc['BP']}
        objgoeabp = _get_objgoeans(pop, bp2assc, godag)
        results_bp = objgoeabp.run_study(enriched)
        print('{N} results'.format(N=len(results_bp)))
    print("TEST PASSED")


def _get_objgoeans(pop, ns2assoc, godag):
    """Run gene ontology enrichment analysis (GOEA)."""
    return GOEnrichmentStudyNS(pop, ns2assoc, godag,
                               propagate_counts=True,
                               relationships=False,
                               alpha=0.05,
                               methods={'fdr_bh'})

def _get_enriched_goids(top, godag):
    """Get a set of GO IDs related to specified top term"""
    gosubdag = GoSubDag(None, godag, relationships=True)
    return {go for go, s in gosubdag.rcntobj.go2descendants.items() if top in s or top == go}

def _get_objanno(fin_anno, anno_type=None, **kws):
    """Get association object"""
    full_anno = os.path.join(REPO, fin_anno)
    dnld_annofile(full_anno, anno_type)
    obj = get_objanno(full_anno, anno_type, **kws)
    return obj


if __name__ == '__main__':
    test_find_enrichment()

# Copyright (C) 2010-2019, DV Klopfenstein, H Tang. All rights reserved.
