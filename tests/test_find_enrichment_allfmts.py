#!/usr/bin/env python3
"""Test running an enrichment using any annotation file format."""

from __future__ import print_function

__copyright__ = "Copyright (C) 2010-2019, DV Klopfenstein, H Tang. All rights reserved."

import os
from goatools.base import get_godag
# from goatools.cli.find_enrichment import GoeaCliFnc
# from goatools.test_data.cli.find_enrichment_dflts import ArgsDict
from goatools.associations import dnld_annofile
from goatools.anno.factory import get_objanno
from goatools.gosubdag.gosubdag import GoSubDag
from goatools.go_enrichment import GOEnrichmentStudy

REPO = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..")


def test_find_enrichment():
    """RUn an enrichments using all annotation file formats"""

    godag = get_godag("go-basic.obo", optional_attrs=['relationship'])
    gos = _get_enriched_goids('GO:0006959', godag)  # GO IDs related to humoral response

    # pylint: disable=superfluous-parens
    print('- DOWNLOAD AND LOAD -----------------------------------------------')
    annoobjs = [
        _get_objanno('gene2go', taxid=10090),
        _get_objanno('goa_human.gaf'),
        _get_objanno('goa_human.gpad'),
        _get_objanno('data/association', anno_type='id2gos'),
    ]

    cmdline_args = []
    for obj in annoobjs:
        assc = obj.get_id2gos()
        pop = obj.get_population()
        enriched = obj.get_ids_g_goids(gos)
        fout_pop = 'ids_pop_{BASE}.txt'.format(BASE=obj.name)
        fout_stu = 'ids_stu_{BASE}.txt'.format(BASE=obj.name)
        _wr(fout_pop, pop)
        _wr(fout_stu, list(enriched)[:100])
        cmdline_args.append((fout_stu, fout_pop, obj.filename))
        objgoea = _get_objgoea(pop, assc, godag)
        results = objgoea.run_study(enriched)
    print("TEST PASSED")

    for args in cmdline_args:
        print('scripts/find_enrichment.py {FILES}'.format(FILES=' '.join(args)))


def _wr(fout_txt, genes):
    """Write genes into a text file."""
    with open(fout_txt, 'w') as prt:
        for gene in sorted(genes):
            prt.write('{GENE}\n'.format(GENE=gene))
        print('  {N:6,} WROTE: {TXT}'.format(N=len(genes), TXT=fout_txt))

def _get_objgoea(pop, assoc, godag):
    """Run gene ontology enrichment analysis (GOEA)."""
    return GOEnrichmentStudy(pop, assoc, godag,
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
