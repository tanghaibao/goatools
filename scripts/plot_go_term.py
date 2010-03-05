#!/usr/bin/env python
# -*- coding: UTF-8 -*-

import os
import sys
from goatools.obo_parser import GODag


if __name__ == '__main__':

    import optparse
    p = optparse.OptionParser("%prog [obo_file]")
    p.add_option("--description", dest="desc", 
            help="write term descriptions to stdout" \
                 " from the obo file specified in args", action="store_true")
    p.add_option("--term", dest="term", help="write the parents and children" \
            "of the query term", action="store", type="string", default=None)

    (options, args) = p.parse_args()

    if not len(args):
        obo_file = None
    else:
        obo_file = args[0]
        assert os.path.exists(obo_file), "file %s not found!" % obo_file

    if obo_file is None:
        g = GODag()
    else:
        g = GODag(obo_file)

    if options.desc:
        g.write_dag()

    # run a test case
    if options.term is not None:
        g.query_term(options.term, draw_lineage=True, verbose=True)

