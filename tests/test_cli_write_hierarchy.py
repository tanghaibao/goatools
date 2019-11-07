#!/usr/bin/env python
"""Test that hierarchy below specified GO terms is printed."""

from __future__ import print_function

__copyright__ = "Copyright (c) 2017-present, DV Klopfenstein. Haiboa Tang. All rights reserved."

import os
from goatools.cli.wr_hierarchy import WrHierCli
from goatools.base import download_go_basic_obo
REPO = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..")

#  --o               Output file in ASCII text format
#  --no_indent       Do not indent GO terms
#  --max_indent      max depth for printing relative to GO Term
#  --num_child       Print count of total number of children for each GO
#  --concise         If a branch has already been printed, do not re-print.
#                    Print '===' instead of dashes to note the point of compression

def test_cli():
    """Add and remove markers for a file."""
    # pylint: disable=bad-whitespace
    args_exp = [
        # args                   exp_set expected_dict
        # --------               ------- ---------------------
        ([],                     {'dag':'go-basic.obo', 'dash_len':6}),
        (['--dag=go-basic.obo'], {'dag':'go-basic.obo', 'dash_len':6}),
        (['-o rpt.txt'],         {'dag':'go-basic.obo', 'dash_len':6, 'o':'rpt.txt'}),
        (['--max_indent=7'],     {'dag':'go-basic.obo', 'dash_len':6, 'max_indent':7}),
        (['CC', '--concise'],    {'dag':'go-basic.obo', 'dash_len':6, 'GO':['CC'], 'concise':True}),
        (['--no_indent'],        {'dag':'go-basic.obo', 'dash_len':6, 'no_indent':True}),
        (['--concise', '--no_indent'], {'dag':'go-basic.obo', 'dash_len':6,
                                        'concise':True, 'no_indent':True}),
    ]
    download_go_basic_obo('go-basic.obo', loading_bar=None)
    for idx, (args, exp_dict) in enumerate(args_exp):
        print("ARGS={ARGS}".format(ARGS=args))
        print("EXP={EXP}".format(EXP=exp_dict))
        obj = WrHierCli(args)
        print("DCT: {DCT}".format(DCT=obj.kws))
        print("WWWWWWWWWWWWWWWWWWW WrHierCli", obj.kws)
        assert obj.kws == exp_dict, "DCT: ACT({}) != EXP({})".format(obj.kws, exp_dict)
        print("")
        # Test writing to a file
        if obj.goids:
            fout_txt = os.path.join(REPO, 'wrhier{N}.txt'.format(N=idx))
            os.system('rm -f {FILE}'.format(FILE=fout_txt))
            obj.wrtxt_hier(fout_txt)
            assert os.path.exists(fout_txt), 'FILE NOT FOUND({F})'.format(F=fout_txt)


if __name__ == '__main__':
    test_cli()

# Copyright (c) 2017-present, DV Klopfenstein, Haibao Tang. All rights reserved.
