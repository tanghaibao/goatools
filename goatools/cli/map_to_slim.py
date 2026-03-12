"""Map GO terms or associations to GO slim terms."""

from __future__ import print_function

__copyright__ = "Copyright (C) 2010-present, H Tang et al. All rights reserved."
__author__ = "various"

import argparse
import os
import sys

from ..associations import read_associations
from ..mapslim import mapslim
from ..obo_parser import GODag


def _get_argparser():
    """Return the CLI argument parser."""
    parser = argparse.ArgumentParser(
        prog="map_to_slim",
        description=__doc__,
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument("go_obo_file", help="Path to the source GO OBO file")
    parser.add_argument("goslim_obo_file", help="Path to the GO slim OBO file")
    parser.add_argument(
        "--term",
        dest="term",
        help="GO term to map to GO slim IDs. Mutually exclusive with --association_file",
    )
    parser.add_argument(
        "--association_file",
        dest="ass_file_name",
        help=(
            "Protein products and their associations to map to GO slim terms. "
            "Mutually exclusive with --term"
        ),
    )
    parser.add_argument(
        "--slim_out",
        dest="slim_out",
        default="direct",
        choices=("direct", "all"),
        help="Report direct slim ancestors or all slim ancestors",
    )
    return parser


def _read_dags(args):
    """Read the GO DAG and GO slim DAG."""
    for fin in (args.go_obo_file, args.goslim_obo_file):
        assert os.path.exists(fin), "file %s not found!" % fin
    return GODag(args.go_obo_file), GODag(args.goslim_obo_file)


def _write_term(go_term, go_dag, goslim_dag, only_direct, prt):
    """Write GO slim mappings for a single GO term."""
    if go_term not in go_dag:
        sys.stderr.write("term %s not found!\n" % go_term)
        raise SystemExit(1)
    direct_anc, all_anc = mapslim(go_term, go_dag, goslim_dag)
    slim_terms = direct_anc if only_direct else all_anc
    prt.write("{SLIMS}\n".format(SLIMS=";".join(sorted(slim_terms))))


def _write_associations(fin_assc, go_dag, goslim_dag, only_direct, prt):
    """Write GO slim mappings for all gene associations in a file."""
    assert os.path.exists(fin_assc), "file %s not found!" % fin_assc
    assocs = read_associations(fin_assc, "id2gos")
    for protein_product, go_terms in assocs.items():
        all_covered_anc = set()
        all_all_anc = set()
        for go_term in go_terms:
            if go_term not in go_dag:
                continue
            direct_anc, all_anc = mapslim(go_term, go_dag, goslim_dag)
            all_all_anc |= all_anc
            all_covered_anc |= all_anc - direct_anc
        slim_terms = all_all_anc - all_covered_anc if only_direct else all_all_anc
        prt.write(
            "{PROD}\t{SLIMS}\n".format(
                PROD=protein_product, SLIMS=";".join(sorted(slim_terms))
            )
        )


def main(args=None):
    """Run the GO slim mapping CLI."""
    parser = _get_argparser()
    nspc = parser.parse_args(args)
    if (nspc.term is None) == (nspc.ass_file_name is None):
        parser.error("Specify exactly one of --term or --association_file")
    go_dag, goslim_dag = _read_dags(nspc)
    only_direct = nspc.slim_out == "direct"
    if nspc.term:
        _write_term(nspc.term, go_dag, goslim_dag, only_direct, sys.stdout)
        return
    _write_associations(nspc.ass_file_name, go_dag, goslim_dag, only_direct, sys.stdout)
