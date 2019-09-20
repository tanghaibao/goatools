from __future__ import absolute_import

#### # make the module mportable
#### from goatools.go_enrichment import *
#### from . import multiple_testing
#### from . import obo_parser

__author__ = ("Haibao Tang", "DV Klopfenstein")
__copyright__ = "Copyright (C) 2009-present, Haibao Tang, DV Klopfenstein"
__email__ = "tanghaibao@gmail.com"
__license__ = "BSD"
__status__ = "Development"

from ._version import get_versions
__version__ = get_versions()['version']
del get_versions
