#!/usr/bin/env python
# -*- coding: UTF-8 -*-

import os
import os.path as op
import sys
sys.path.insert(0, op.join(op.dirname(__file__), ".."))
from goatools.obo_parser import GODag


if __name__ == '__main__':

    import optparse
    p = optparse.OptionParser("%prog [go_obo_file] [goslim_obo_file]")
    p.add_option("--term", dest="term", help="a term (accession id) to map "
                 "to slim IDs", action="store", type="string", default=None)
    p.add_option("--accession_file", dest="ass_file_name", action="store",
                 help="TODO TODO TODO", type="string", default=None)

    opts, args = p.parse_args()

    usage = ("Usage: map_to_slim --term go_term "
             "<go_obo_file> <goslim_obo_file>\n"
             "or   : map_to_slim --accession_file filename "
             "<go_obo_file> <goslim_obo_file>\n")

    # check for correct number of arguments
    if len(args) != 2:
        print >>sys.stderr, usage
        sys.exit(1)

    obo_file = args[0]
    assert os.path.exists(obo_file), "file %s not found!" % obo_file

    slim_obo_file = args[1]
    assert os.path.exists(slim_obo_file), "file %s not found!" % slim_obo_file

    # check that either --term or --accession_file is set
    if (opts.term is None and opts.ass_file_name is None) \
            or ((not opts.term is None) and (not opts.ass_file_name is None)):
        print >>sys.stderr, usage
        sys.exit(1)

    # load DAGs
    go_dag = GODag(obo_file)
    goslim_dag = GODag(slim_obo_file)

    # TODO: read term or accessions from file
    #       and execute the mapslim algorithm per input
    pass
