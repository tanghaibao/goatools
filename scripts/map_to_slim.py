#!/usr/bin/env python
# -*- coding: UTF-8 -*-

import sys

from goatools.cli.main import main as goatools_main


if __name__ == '__main__':
    goatools_main(["map_to_slim"] + sys.argv[1:])
