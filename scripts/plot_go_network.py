#!/usr/bin/env python
# -*- coding: UTF-8 -*-

import os
import sys
from goatools import GODag

if __name__ == '__main__':

    import optparse
    p = optparse.OptionParser("%prog infile")
    (options, args) = p.parse_args()

    if not len(args):
        sys.exit(p.print_help())

    infile = args[0]
    assert os.path.exists(infile), "file %s not found!" % infile

    g = GODag()

