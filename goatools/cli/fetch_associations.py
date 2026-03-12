"""Download associations from GOlr for a given species."""

from __future__ import print_function

__copyright__ = "Copyright (C) 2010-present, H Tang et al. All rights reserved."
__author__ = "various"

import argparse
import sys


def _get_argparser():
    """Return the CLI argument parser."""
    parser = argparse.ArgumentParser(
        prog="fetch_associations",
        description=__doc__,
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument(
        "--taxon_id",
        required=True,
        help=(
            "NCBI taxon ID. This must match the exact species or strain used by GO Central"
        ),
    )
    parser.add_argument(
        "--golr_url",
        default="http://golr.geneontology.org/solr/",
        help="GOlr endpoint URL",
    )
    parser.add_argument(
        "-o",
        default=None,
        help="Optional output filename. Defaults to stdout",
    )
    parser.add_argument(
        "--max_rows",
        default=100000,
        type=int,
        help="Maximum rows to fetch",
    )
    return parser


def main(args=None):
    """Fetch simple gene-term associations from GOlr."""
    try:
        import pysolr
    except ImportError as exc:
        raise SystemExit(
            "fetch_associations requires the optional dependency 'pysolr'"
        ) from exc

    nspc = _get_argparser().parse_args(args)
    solr = pysolr.Solr(nspc.golr_url, timeout=30)

    sys.stderr.write("TAX:{TAX}\n".format(TAX=nspc.taxon_id))
    results = solr.search(
        q='document_category:"bioentity" AND taxon:"NCBITaxon:{TAX}"'.format(
            TAX=nspc.taxon_id
        ),
        fl="bioentity_label,annotation_class_list",
        rows=nspc.max_rows,
    )
    sys.stderr.write("NUM GENES:{N}\n".format(N=len(results)))
    if len(results) == 0:
        sys.stderr.write("NO RESULTS")
        raise SystemExit(1)
    if len(results) == nspc.max_rows:
        sys.stderr.write("max_rows set too low")
        raise SystemExit(1)

    file_out = sys.stdout if nspc.o is None else open(nspc.o, "w")
    try:
        for result in results:
            gene_symbol = result["bioentity_label"]
            sys.stderr.write("{GENE}\n".format(GENE=gene_symbol))
            if "annotation_class_list" in result:
                file_out.write(
                    "{GENE}\t{GOS}\n".format(
                        GENE=gene_symbol,
                        GOS=";".join(result["annotation_class_list"]),
                    )
                )
            else:
                sys.stderr.write("no annotations for {GENE}\n".format(GENE=gene_symbol))
    finally:
        if nspc.o is not None:
            file_out.close()

    if nspc.o is not None:
        sys.stdout.write("  WROTE: {}\n".format(nspc.o))
