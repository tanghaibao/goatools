"""Used to find all genes or gene products annotated w/GO terms that match a regex."""

import sys
from goatools.obo_parser import GODag

__copyright__ = "Copyright (C) 2010-2018, DV Klopfenstein, H Tang, All rights reserved."
__author__ = "DV Klopfenstein"

class GoSearch(object):
    """Returns GOs matching a regex pattern."""

    def __init__(self, fin_go_basic_obo, go2items, log=None):
        self.log = sys.stdout if log is None else log
        self.bstdout = True if log is None else log
        # Some obo fields often used in searching. Many are optional to load when reading obo
        self.goa_srch_hdrs = ['defn', 'comment', 'name', 'is_a', 'relationship', 'synonym', 'xref']
        self.obo_dag = GODag(fin_go_basic_obo, optional_attrs=self.goa_srch_hdrs)
        self.go2items = go2items

    def get_matching_gos(self, compiled_pattern, **kws):
        """Return all GOs which match the user regex pattern."""
        # kws: prt gos
        matching_gos = []
        obo_dag = self.obo_dag
        prt = kws['prt'] if 'prt' in kws else self.log
        prt.write('\nPATTERN SEARCH: "{P}"\n'.format(P=compiled_pattern.pattern))
        # Only look through GOs in annotation or user-specified GOs
        srchgos = kws['gos'] if 'gos' in kws else self.go2items.keys()
        for go_id in srchgos:
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
        # Print summary message
        self._summary_matching_gos(prt, compiled_pattern.pattern, matching_gos, srchgos)
        return matching_gos

    @staticmethod
    def _summary_matching_gos(prt, pattern, matching_gos, all_gos):
        """Print summary for get_matching_gos."""
        msg = 'Found {N} GO(s) out of {M} matching pattern("{P}")\n'
        num_gos = len(matching_gos)
        num_all = len(all_gos)
        prt.write(msg.format(N=num_gos, M=num_all, P=pattern))

    def _search_vals(self, compiled_pattern, fld_val):
        """Search for user-regex in scalar or iterable data values."""
        matches = []
        if isinstance(fld_val, set):
            for val in fld_val:
                self._search_val(matches, compiled_pattern, val)
        elif isinstance(fld_val, str):
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

# Copyright (C) 2010-2018, DV Klopfenstein, H Tang, All rights reserved.
