#!/usr/bin/env python
# -*- coding: UTF-8 -*-

from __future__ import print_function
import os
import os.path as op
import sys
sys.path.insert(0, op.join(op.dirname(__file__), ".."))
from goatools.obo_parser import GODag
from goatools.mapslim import mapslim
from goatools.associations import read_associations


if __name__ == '__main__':

    import optparse
    p = optparse.OptionParser("%prog [options] go_obo_file goslim_obo_file")
    p.add_option("--term", dest="term", help="a term (association id) to map "
                 "to slim IDs. This can not be used together with "
                 "--association_file", action="store", type="string",
                 default=None)
    p.add_option("--association_file", dest="ass_file_name", action="store",
                 help="the file of protein products and their associations "
                 "to be mapped to GO slim terms. This can not be used "
                 "together with --term", type="string", default=None)
    p.add_option("--slim_out", dest="slim_out", action="store", type="string",
                 default="direct", help="One of `direct` or `all`. Defines "
                 "whether the output should contain all slim terms (all "
                 "ancestors) or only direct slim terms (only direct "
                 "ancestors)")
    opts, args = p.parse_args()

    # check for correct number of arguments
    if len(args) != 2:
        p.print_help()
        sys.exit(1)

    obo_file = args[0]
    assert os.path.exists(obo_file), "file %s not found!" % obo_file

    slim_obo_file = args[1]
    assert os.path.exists(slim_obo_file), "file %s not found!" % slim_obo_file

    # check that either --term or --association_file is set
    if (opts.term is None and opts.ass_file_name is None) \
            or ((opts.term is not None) and (opts.ass_file_name is not None)):
        p.print_help()
        sys.exit(1)

    # check that slim_out is either "direct" or "all" and set according flag
    only_direct = None
    if opts.slim_out == "direct":
        only_direct = True
    elif opts.slim_out == "all":
        only_direct = False
    else:
        p.print_help()
        sys.exit(1)

    # load DAGs
    go_dag = GODag(obo_file)
    goslim_dag = GODag(slim_obo_file)

    # in case a single term is given as input:
    if opts.term:
        if opts.term not in go_dag:
            print(("term %s not found!" % opts.term), file=sys.stderr)
            sys.exit(1)
        direct_anc, all_anc = mapslim(opts.term, go_dag, goslim_dag)
        # output either all or only direct slims, depending on user command
        if only_direct:
            slim_terms_str = ";".join(direct_anc)
        else:
            slim_terms_str = ";".join(all_anc)
        print(slim_terms_str)

    # in case a association file is given as input
    if opts.ass_file_name:
        assert os.path.exists(opts.ass_file_name), ("file %s not found!"
                                                    % opts.ass_file_name)
        assocs = read_associations(opts.ass_file_name)
        for protein_product, go_terms in assocs.items():
            all_direct_anc = set()
            all_covered_anc = set()
            all_all_anc = set()
            for go_term in go_terms:
                if go_term not in go_dag:
                    continue
                direct_anc, all_anc = mapslim(go_term, go_dag, goslim_dag)
                all_all_anc |= all_anc
                # collect all covered ancestors, so the direct ancestors
                # can be calculated afterwards
                all_covered_anc |= (all_anc - direct_anc)
            all_direct_anc = all_all_anc - all_covered_anc
            # output either all or only direct, depending on user command
            if only_direct:
                slim_terms_str = ";".join(all_direct_anc)
            else:
                slim_terms_str = ";".join(all_all_anc)
            print((protein_product + "\t" + slim_terms_str))
