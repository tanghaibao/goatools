#!/usr/bin/env python
# -*- coding: UTF-8 -*-
from __future__ import print_function
import sys
try:
    from exceptions import EOFError
except ImportError:
    pass

from goatools.obo_parser import GODag

class GODagInfo:

    def __init__(self, obo_file="go-basic.obo"):

        self.obo_dag = GODag(obo_file)
        self.nslst = ['biological_process', 'molecular_function', 'cellular_component']
        self.max_vals = self._mk_max_vals()
        self.prt_maxvals(sys.stdout)

    def prt_maxvals(self, out=sys.stdout):
        """Prints max values of both level and depth."""
        for desc, nss in self.max_vals.items():
          for ns in self.nslst:
            out.write('{} {} {}\n'.format(desc, ns, nss[ns]))
 
    def _mk_max_vals(self):
        """Initialize the max value for BP, MF, CC for both depth and level.

           Example max depth and level values for BP, MF and CC:

               biological_process depth 16
               molecular_function depth 15
               cellular_component depth 12

               biological_process level 12
               molecular_function level 11
               cellular_component level 9
        """
        max_vals = {}
        num_namespaces = len(self.nslst)
        cnts = self.obo_dag.get_cnts_levels_depths_recs(set(self.obo_dag.values()))
        for desc, cnts_curr in cnts.items():
          max_vals[desc] = {}
          for lev, ctr in sorted(cnts_curr.items(), reverse=True):
            if len(max_vals[desc]) != num_namespaces:
              for ns, val in ctr.items():
                if ns not in max_vals[desc]:
                  max_vals[desc][ns] = lev
            else:
              break
        return max_vals
  


