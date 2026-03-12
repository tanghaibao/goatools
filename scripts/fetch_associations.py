"""Script for downloading associations from GO Database for a given species."""
# -*- coding: UTF-8 -*-

import sys

from goatools.cli.main import main as goatools_main

if __name__ == "__main__":
    goatools_main(["fetch_associations"] + sys.argv[1:])
