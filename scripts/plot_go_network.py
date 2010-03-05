#!/usr/bin/env python
# -*- coding: UTF-8 -*-

import os
import sys
from goatools import GODag


def read_go_counts(infile):

    fp = file(infile)
    term_cnt = {}
    for row in fp:
        term, cnt = row.split()
        term_cnt[term] = int(cnt)
    return term_cnt


if __name__ == '__main__':

    import optparse
    p = optparse.OptionParser("%prog infile")
    (options, args) = p.parse_args()

    if not len(args):
        sys.exit(p.print_help())

    infile = args[0]
    assert os.path.exists(infile), "file %s not found!" % infile

    g = GODag()
    term_cnt = read_go_counts(infile)

