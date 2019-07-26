#!/usr/bin/env python3
"""Test running an enrichment using any annotation file format."""

from __future__ import print_function

__copyright__ = "Copyright (C) 2010-2019, DV Klopfenstein, H Tang. All rights reserved."

import os
import sys

# REPO = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..")


def test_find_enrichment(run_all=False):
    """RUn an enrichments using all annotation file formats"""

    if run_all:
        for idx, cmd in enumerate(_get_cmds()):
            print('------------------- TEST {I} ------------------------------------'.format(I=idx))
            print('CMD: {CMD}'.format(CMD=cmd))
            assert os.system(cmd) == 0
        print("TEST PASSED")
    else:
        print('RUN THIS TEST WITH AN ARGUMENT')


def _get_cmds():
    """Get commands used in ./doc/md/README_find_enrichment.md"""
    # pylint: disable=line-too-long
    return [
        'python3 scripts/find_enrichment.py data/study data/population data/association --outfile=goea.xlsx,goea.tsv --pval_field=fdr_bh',
        'python3 scripts/find_enrichment.py data/study data/population data/association --pval=0.05 --method=fdr_bh --pval_field=fdr_bh --outfile=results_id2gos.xlsx',
        'python3 scripts/find_enrichment.py ids_stu_gaf.txt ids_pop_gaf.txt goa_human.gaf --pval=0.05 --method=fdr_bh --pval_field=fdr_bh --outfile=results_gaf.xlsx',
        'python3 scripts/find_enrichment.py ids_stu_gpad.txt ids_pop_gpad.txt goa_human.gpad --pval=0.05 --method=fdr_bh --pval_field=fdr_bh --outfile=results_gpad.xlsx',
        'python3 scripts/find_enrichment.py ids_stu_gpad.txt ids_pop_gpad.txt goa_human.gpad --pval=0.05 --method=fdr_bh --pval_field=fdr_bh --outfile=results_gpad.xlsx -r',
        'python3 scripts/find_enrichment.py ids_stu_gpad.txt ids_pop_gpad.txt goa_human.gpad --pval=0.05 --method=fdr_bh --pval_field=fdr_bh --outfile=results_gpad.xlsx --relationships=regulates,negatively_regulates,positively_regulates',
        'python3 scripts/find_enrichment.py ids_stu_gene2go_9606.txt ids_pop_gene2go_9606.txt gene2go --pval=0.05 --method=fdr_bh --pval_field=fdr_bh --outfile=results_gene2go_9606.xlsx',
        'python3 scripts/find_enrichment.py ids_stu_gene2go_10090.txt ids_pop_gene2go_10090.txt gene2go --taxid=10090 --pval=0.05 --method=fdr_bh --pval_field=fdr_bh --outfile=results_gene2go_10090.xlsx',
        'python3 scripts/find_enrichment.py ids_stu_gpad.txt ids_pop_gpad.txt goa_human.gpad --ev_exc=IEA --pval=0.05 --method=fdr_bh --pval_field=fdr_bh --outfile=results_gpad.xlsx',
        'python3 scripts/find_enrichment.py ids_stu_gpad.txt ids_pop_gpad.txt goa_human.gpad --ev_inc=Experimental --pval=0.05 --method=fdr_bh --pval_field=fdr_bh --outfile=results_gpad.xlsx',
        'python3 scripts/find_enrichment.py ids_stu_gpad.txt ids_pop_gpad.txt goa_human.gpad --ev_inc=EXP,IDA,IPI,IMP,IGI,IEP --pval=0.05 --method=fdr_bh --pval_field=fdr_bh --outfile=results_gpad.xlsx',
        'python3 scripts/find_enrichment.py --ev_help',
        'python3 scripts/find_enrichment.py data/study data/population data/association',
        'python3 scripts/find_enrichment.py data/study data/population data/association --sections=goatools.test_data.sections.data2018_07_find_enrichment',
        'python3 scripts/find_enrichment.py data/study data/population data/association --outfile=goea_uncorr.xlsx',
        'python3 scripts/find_enrichment.py data/study data/population data/association --outfile=goea_fdr_bh.xlsx,goea_fdr_bh.tsv --pval_field=fdr_bh',
        'python3 scripts/find_enrichment.py data/study data/population data/association --outfile=goea_fdr_bh_flat.xlsx --method=fdr_bh',
        'python3 scripts/find_enrichment.py data/study data/population data/association --outfile=goea_fdr_bh_grpd.xlsx --method=fdr_bh --sections=goatools.test_data.sections.data2018_07_find_enrichment',
        'python3 scripts/find_enrichment.py data/study data/population data/association --outfile=goea_all.xlsx,goea_all.tsv --pval=-1',
    ]


if __name__ == '__main__':
    test_find_enrichment(len(sys.argv) != 1)

# Copyright (C) 2010-2019, DV Klopfenstein, H Tang. All rights reserved.
