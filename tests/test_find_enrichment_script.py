#!/usr/bin/env python3
"""Test running an enrichment using any annotation file format."""

from __future__ import print_function

__copyright__ = "Copyright (C) 2010-2019, DV Klopfenstein, H Tang. All rights reserved."

import os
from itertools import chain
import collections as cx
from goatools.base import get_godag
from goatools.associations import dnld_annofile
from goatools.anno.factory import get_objanno
from goatools.gosubdag.gosubdag import GoSubDag

REPO = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..")


# pylint: disable=fixme,too-many-locals
def test_find_enrichment():
    """RUn an enrichments using all annotation file formats"""

    godag = get_godag("go-basic.obo", optional_attrs=['relationship'])
    e_goids = _get_enriched_e_goids('GO:0006959', godag)  # GO IDs related to humoral response

    # pylint: disable=superfluous-parens
    print('- DOWNLOAD AND LOAD -----------------------------------------------')
    annoobjs = [
        _get_objanno('gene2go', taxid=10090),
        _get_objanno('gene2go', taxid=9606),
        _get_objanno('goa_human.gaf'),
        _get_objanno('goa_human.gpad', godag=godag),
        _get_objanno('data/association', anno_type='id2gos', godag=godag),
    ]

    pat = ('python3 scripts/find_enrichment.py {STU} {POP} {ASSC} '
           '--pval=0.05 --method=fdr_bh --pval_field=fdr_bh '
           '--taxid={TAXID} {INC} {EXC} --outfile=results_{NAME}.xlsx')
    cmds = []
    for obj in annoobjs:
        ns2assc = obj.get_ns2assc()
        _idngos_list = list(chain.from_iterable([k2v.items() for k2v in ns2assc.values()]))
        pop = set(d for d, _ in _idngos_list)
        # TODO: 20,263 pop IDs      6,847 stu IDs      2,884 int IDs
        enriched = set(nt.DB_ID for nt in obj.get_associations() if nt.GO_ID in e_goids)
        stu = enriched.intersection(pop)
        print('{N:6,} pop IDs: {ID}'.format(N=len(pop), ID=list(pop)[:4]))
        print('{N:6,} enr IDs: {ID}'.format(N=len(enriched), ID=list(enriched)[:4]))
        print('{N:6,} int IDs: {ID}'.format(N=len(stu), ID=list(stu)[:4]))
        fout_pop = os.path.join(REPO, 'ids_pop_{BASE}.txt'.format(BASE=obj.get_name()))
        fout_stu = os.path.join(REPO, 'ids_stu_{BASE}.txt'.format(BASE=obj.get_name()))
        _wr(fout_pop, pop)
        _wr(fout_stu, list(stu)[:100])
        cmd = pat.format(STU=fout_stu, POP=fout_pop, ASSC=obj.filename,
                         TAXID=obj.get_taxid(), NAME=obj.get_name(),
                         INC='', EXC='')
        cmds.append(cmd)
        print('\nRUNNING {NAME}: {CMD}\n'.format(CMD=cmd, NAME=obj.get_name()))
        assert os.system(cmd) == 0

    fout_scr = 'test_find_enrichment_script.sh'
    with open(fout_scr, 'w') as prt:
        print("COMANDS RUN:")
        for cmd in cmds:
            print(cmd)
            prt.write('{CMD}\n'.format(CMD=cmd))
        print('  WROTE: {SCRIPT}'.format(SCRIPT=fout_scr))


    print("TEST PASSED")


def _wr(fout_txt, genes):
    """Write genes into a text file."""
    with open(fout_txt, 'w') as prt:
        for gene in sorted(genes):
            prt.write('{GENE}\n'.format(GENE=gene))
        print('  {N:6,} WROTE: {TXT}'.format(N=len(genes), TXT=fout_txt))

def _get_enriched_e_goids(top, godag):
    """Get a set of GO IDs related to specified top term"""
    gosubdag = GoSubDag(None, godag, relationships=True)
    e_goids = {go for go, s in gosubdag.rcntobj.go2descendants.items() if top in s or top == go}
    ctr = cx.Counter([gosubdag.go2nt[go].NS for go in e_goids])
    print('{N:6,} Enriched GO IDs: {CTR}'.format(N=len(e_goids), CTR=ctr.most_common()))
    return e_goids

def _get_objanno(fin_anno, anno_type=None, **kws):
    """Get association object"""
    full_anno = os.path.join(REPO, fin_anno)
    dnld_annofile(full_anno, anno_type)
    obj = get_objanno(full_anno, anno_type, **kws)
    return obj


if __name__ == '__main__':
    test_find_enrichment()

# Copyright (C) 2010-2019, DV Klopfenstein, H Tang. All rights reserved.
