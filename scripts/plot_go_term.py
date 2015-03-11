#!/usr/bin/env python
# -*- coding: UTF-8 -*-

import os
import os.path as op
import sys
sys.path.insert(0, op.join(op.dirname(__file__), ".."))
from goatools.obo_parser import GODag, GraphEngines


if __name__ == '__main__':

    import optparse
    p = optparse.OptionParser("%prog [obo_file]")
    p.add_option("--description", dest="desc",
                 help="write term descriptions to stdout"
                 " from the obo file specified in args", action="store_true")
    p.add_option("--term", dest="term", help="write the parents and children"
                 " of the query term", action="store", type="string",
                 default=None)
    p.add_option("--engine", default="pygraphviz",
                 choices=GraphEngines,
                 help="Graph plot engine, must be one of {0} [default: %default]".\
                        format("|".join(GraphEngines))
                 )
    p.add_option("--gml", action="store_true",
                 help="Write GML output (for Cytoscape) [default: %default]")
    p.add_option("--disable-draw-parents",
                 action="store_false",
                 dest='draw_parents',
                 help="Do not draw parents of the query term")
    p.add_option("--disable-draw-children",
                 action="store_false",
                 dest='draw_children',
                 help="Do not draw children of the query term")

    p.set_defaults(draw_parents=True)
    p.set_defaults(draw_children=True)

    opts, args = p.parse_args()

    if not len(args):
        obo_file = "go-basic.obo"
    else:
        obo_file = args[0]
        assert os.path.exists(obo_file), "file %s not found!" % obo_file

    g = GODag(obo_file)

    if opts.desc:
        g.write_dag()

    # run a test case
    if opts.term is not None:
        rec = g.query_term(opts.term, verbose=True)
        g.draw_lineage([rec], engine=opts.engine,
                       gml=opts.gml,
                       draw_parents=opts.draw_parents,
                       draw_children=opts.draw_children)
