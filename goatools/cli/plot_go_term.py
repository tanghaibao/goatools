"""Plot a GO term lineage using the GODag drawing helpers."""

__copyright__ = "Copyright (C) 2010-present, H Tang et al. All rights reserved."
__author__ = "various"

import argparse
import os

from ..obo_parser import GODag, GraphEngines


def _get_argparser():
    """Return the CLI argument parser."""
    parser = argparse.ArgumentParser(
        prog="goatools plot_go_term",
        description=__doc__,
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument(
        "obo_file",
        nargs="?",
        default="go-basic.obo",
        help="Path to the GO OBO file",
    )
    parser.add_argument(
        "--description",
        dest="desc",
        action="store_true",
        help="Write term descriptions to stdout from the specified OBO file",
    )
    parser.add_argument(
        "--term",
        dest="term",
        help="Write the parents and children of the query term",
    )
    parser.add_argument(
        "--engine",
        default="pygraphviz",
        choices=GraphEngines,
        help="Graph plot engine",
    )
    parser.add_argument(
        "--gml",
        action="store_true",
        help="Write GML output for Cytoscape",
    )
    parser.add_argument(
        "--disable-draw-parents",
        action="store_false",
        dest="draw_parents",
        help="Do not draw parents of the query term",
    )
    parser.add_argument(
        "--disable-draw-children",
        action="store_false",
        dest="draw_children",
        help="Do not draw children of the query term",
    )
    parser.add_argument(
        "--output",
        "-o",
        default="GO_lineage.pdf",
        help="Output filename. The suffix determines the image format",
    )
    parser.add_argument(
        "--dpi",
        default=96,
        type=int,
        help="Figure DPI for raster outputs",
    )
    parser.add_argument(
        "--wrap-width",
        default=999,
        type=int,
        help="Maximum width of graph nodes in characters",
    )
    parser.set_defaults(draw_parents=True, draw_children=True)
    return parser


def main(args=None):
    """Run the GO term lineage plotting CLI."""
    parser = _get_argparser()
    opts = parser.parse_args(args)
    if not os.path.exists(opts.obo_file):
        parser.error("file %s not found!" % opts.obo_file)
    godag = GODag(opts.obo_file)
    if opts.desc:
        godag.write_dag()
    if opts.term is not None:
        rec = godag.query_term(opts.term, verbose=True)
        godag.draw_lineage(
            [rec],
            dpi=opts.dpi,
            engine=opts.engine,
            gml=opts.gml,
            output=opts.output,
            draw_parents=opts.draw_parents,
            draw_children=opts.draw_children,
            wrap_width=opts.wrap_width,
        )
