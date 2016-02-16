
import sys
import os
from goatools.obo_parser import GODag

"""Used to find all genes or gene products annotated w/GO terms that match a regex."""

ROOT = os.path.dirname(os.path.abspath(__file__)) + "/data/"

__copyright__ = "Copyright (C) 2010-2016, H Tang et al., All rights reserved."
__author__ = "DV Klopfenstein"

class GoSearch(object):
    """Returns GOs matching a regex pattern."""

    def __init__(self, fin_go_basic_obo, go2items, log=sys.stdout):
        self.log = log
        # Some obo fields often used in searching. Many are optional to load when reading obo
        self.goa_srch_hdrs = ['defn', 'comment', 'name', 'is_a', 'relationship', 'synonym', 'xref']
        self.obo_dag = GODag(fin_go_basic_obo, optional_attrs=self.goa_srch_hdrs)
        self.go2items = go2items

    def get_matching_gos(self, compiled_pattern, prt=None):
        """Return all GOs which match the user regex pattern."""
        matching_GOs = []
        obo_dag = self.obo_dag
        if prt is None:
            prt = self.log
        # Only look through GOs in annotation
        for GO in self.go2items:
          oGO = obo_dag.get(GO, None)
          if oGO is not None:
              for hdr in self.goa_srch_hdrs:
                  if hdr in oGO.__dict__:
                      fld_val = getattr(oGO, hdr)
                      matches = self._search_vals(compiled_pattern, fld_val)
                      for m in matches: 
                          prt.write("MATCH {GO}({NAME}) {FLD}: {M}\n".format(
                            FLD=hdr, GO=oGO.id, NAME=oGO.name, M=m))
                      if matches:
                          matching_GOs.append(GO)
          else:
              prt.write("**WARNING: {GO} found in annotation is not found in obo\n".format(GO=GO)) 
        matching_GOs = set(matching_GOs)
        prt.write("{N} GOs out of {M} found for matching pattern({P})\n".format(
            N=len(matching_GOs), M=len(self.go2items), P=compiled_pattern.pattern))
        return matching_GOs

    def _search_vals(self, compiled_pattern, fld_val):
        """Search for user-regex in scalar or iterable data values."""
        matches = []
        if isinstance(fld_val, set):
            for val in fld_val:
                self._search_val(matches, compiled_pattern, val)
        else:
            self._search_val(matches, compiled_pattern, fld_val)
        return matches

    def _search_val(self, matches, compiled_pattern, fld_val):
        """Search for user-regex in scalar data values."""
        M = compiled_pattern.search(fld_val)
        if M:
            matches.append(fld_val)

    def add_children_GOs(self, GOs):
        """Return all GOs who have parents which match the user's compiled regex pattern."""
        lst = []
        obo_dag = self.obo_dag
        rx = lambda oGO: list(oGO.get_all_children()) + [oGO.id]
        for GO in GOs:
            oGO = obo_dag[GO]
            lst.extend(rx(oGO))
        return set(lst)

    def get_items(self, gos):
        """Given GO terms, return genes or gene products for the GOs."""
        items = []
        for go in gos:
            items.extend(self.go2items.get(go, []))
        return set(items)
      
# Copyright (C) 2010-2016, H Tang et al., All rights reserved.
