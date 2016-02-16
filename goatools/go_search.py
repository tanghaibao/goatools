"""Used to find all genes or gene products annotated w/GO terms that match a regex."""

import sys
from goatools.obo_parser import GODag

__copyright__ = "Copyright (C) 2010-2016, DV Klopfenstein, H Tang, All rights reserved."
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
        matching_gos = []
        obo_dag = self.obo_dag
        if prt is None:
            prt = self.log
        # Only look through GOs in annotation
        for go_id in self.go2items:
            go_obj = obo_dag.get(go_id, None)
            if go_obj is not None:
                for hdr in self.goa_srch_hdrs:
                    if hdr in go_obj.__dict__:
                        fld_val = getattr(go_obj, hdr)
                        matches = self._search_vals(compiled_pattern, fld_val)
                        for mtch in matches:
                            prt.write("MATCH {go_id}({NAME}) {FLD}: {M}\n".format(
                                FLD=hdr, go_id=go_obj.id, NAME=go_obj.name, M=mtch))
                        if matches:
                            matching_gos.append(go_id)
            else:
                prt.write("**WARNING: {GO} found in annotation is not found in obo\n".format(
                    GO=go_id))
        matching_gos = set(matching_gos)
        prt.write("{N} GOs out of {M} found for matching pattern({P})\n".format(
            N=len(matching_gos), M=len(self.go2items), P=compiled_pattern.pattern))
        return matching_gos

    def _search_vals(self, compiled_pattern, fld_val):
        """Search for user-regex in scalar or iterable data values."""
        matches = []
        if isinstance(fld_val, set):
            for val in fld_val:
                self._search_val(matches, compiled_pattern, val)
        else:
            self._search_val(matches, compiled_pattern, fld_val)
        return matches

    @staticmethod
    def _search_val(matches, compiled_pattern, fld_val):
        """Search for user-regex in scalar data values."""
        mtch = compiled_pattern.search(fld_val)
        if mtch:
            matches.append(fld_val)

    def add_children_gos(self, gos):
        """Return children of input gos plus input gos."""
        lst = []
        obo_dag = self.obo_dag
        get_children = lambda go_obj: list(go_obj.get_all_children()) + [go_obj.id]
        for go_id in gos:
            go_obj = obo_dag[go_id]
            lst.extend(get_children(go_obj))
        return set(lst)

    def get_items(self, gos):
        """Given GO terms, return genes or gene products for the GOs."""
        items = []
        for go_id in gos:
            items.extend(self.go2items.get(go_id, []))
        return set(items)

# Copyright (C) 2010-2016, DV Klopfenstein, H Tang, All rights reserved.
