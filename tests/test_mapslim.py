#!/usr/bin/env python
# -*- coding: UTF-8 -*-

from __future__ import print_function
import os.path as op
import sys
sys.path.insert(0, op.join(op.dirname(__file__), ".."))
from goatools.obo_parser import GODag
from goatools.mapslim import mapslim

if __name__ == '__main__':
    # read the two mini GO trees (full and slim)
    go_dag = GODag(op.join(op.dirname(__file__), "data/mini_obo.obo"))
    goslim_dag = GODag(op.join(op.dirname(__file__), "data/mini_slim_obo.obo"))

    #
    # This tests the map2slim algorithm with a very small example GO DAG
    # and an even smaller GOslim DAG.
    # The tree and the expected results can be seen at the original
    # map2slim.pl documentation here:
    # http://search.cpan.org/~cmungall/go-perl/scripts/map2slim
    # an image of the graph:
    # http://geneontology.cvs.sourceforge.net/viewvc/geneontology/go-dev/go-perl/doc/map2slim.gif
    #
    # Expected results
    #
    # GO ID  MAPS TO SLIM ID        ALL SLIM ANCESTORS
    # =====  ===============        ==================
    # 5      2+3                    2,3,1
    # 6      3 only                 3,1
    # 7      4 only                 4,3,1
    # 8      3 only                 3,1
    # 9      4 only                 4,3,1
    # 10     2+3                    2,3,1

    expected_results = {
        'GO:0000005': (set(['GO:0000002', 'GO:0000003']),
                       set(['GO:0000001', 'GO:0000002', 'GO:0000003'])),

        'GO:0000006': (set(['GO:0000003']),
                       set(['GO:0000001', 'GO:0000003'])),

        'GO:0000007': (set(['GO:0000004']),
                       set(['GO:0000001', 'GO:0000003', 'GO:0000004'])),

        'GO:0000008': (set(['GO:0000003']),
                       set(['GO:0000001', 'GO:0000003'])),

        'GO:0000009': (set(['GO:0000004']),
                       set(['GO:0000001', 'GO:0000003', 'GO:0000004'])),

        'GO:0000010': (set(['GO:0000002', 'GO:0000003']),
                       set(['GO:0000001', 'GO:0000002', 'GO:0000003']))
    }

    tests_succeed = True

    for go_term, (exp_direct, exp_all) in expected_results.items():
        sys.stderr.write("Testing for term '{}' ...\n".format(go_term))
        direct_anc, all_anc = mapslim(go_term, go_dag, goslim_dag)
        if direct_anc != exp_direct or all_anc != exp_all:
            tests_succeed = False
            sys.stderr.write("failed.\n")
        else:
            sys.stderr.write("success!\n")

    if tests_succeed:
        print("All test passed successfully!")
        sys.exit(0)
    else:
        sys.stderr.write("[ERROR] At least one test failed.\n")
        sys.exit(1)
