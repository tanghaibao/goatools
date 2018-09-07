"""Print a GO term's lower-level hierarchy."""

__copyright__ = "Copyright (C) 2016-2018, DV Klopfenstein, H Tang. All rights reserved."
__author__ = "DV Klopfenstein"

import sys


# pylint: disable=too-many-instance-attributes,too-few-public-methods
class WrHierPrt(object):
    """Print GO hierarchy."""

    def __init__(self, id2obj, id2nt, cfg, prt=sys.stdout):
        self.id2obj = id2obj  # Contains children (and parents)
        self.id2nt = id2nt    # Contains fields for printing
        self.nm2prtfmt = cfg['name2prtfmt']
        self.max_indent = cfg['max_indent']
        self.include_only = cfg['include_only']
        self.go_marks = cfg['go_marks']
        self.concise_prt = cfg['concise_prt']
        self.indent = cfg['indent']
        self.idname = cfg.get('ID', 'id')
        # vars
        self.prt = prt
        self.gos_printed = set()
        self.dash_len = cfg['dash_len'] + 12

    def prt_hier_rec(self, item_id, depth=1):
        """Write hierarchy for a GO Term record and all GO IDs down to the leaf level."""
        ntgo = self.id2nt[item_id]
        ntobj = self.id2obj[item_id]
        # Shortens hierarchy report by only printing the hierarchy
        # for the sub-set of user-specified GO terms which are connected.
        if self.include_only and item_id not in self.include_only:
            return
        nrp = self.concise_prt and item_id in self.gos_printed
        if self.go_marks:
            self.prt.write('{MARK} '.format(MARK='>' if item_id in self.go_marks else ' '))

        # '-' is default character indicating hierarchy level
        # '=' is used to indicate a hierarchical path printed in detail previously.
        dct = ntgo._asdict()
        self.prt.write('{DASHGO:{N}}'.format(
            DASHGO=self._str_dashgoid(dct, depth, not nrp or not ntobj.children),
            N=self.dash_len))

        self.prt.write("{INFO}\n".format(INFO=self.nm2prtfmt['ITEM'].format(**dct)))
        self.gos_printed.add(item_id)
        # Do not print hierarchy below this turn if it has already been printed
        if nrp:
            return
        depth += 1
        if self.max_indent is not None and depth > self.max_indent:
            return
        for child in ntobj.children:
            self.prt_hier_rec(getattr(child, self.idname), depth)

    @staticmethod
    def _str_dash(depth, single_or_double):
        """Return a string containing dashes (optional) and GO ID."""
        # '-' is default character indicating hierarchy level
        # '=' is used to indicate a hierarchical path printed in detail previously.
        letter = '-' if single_or_double else '='
        return ''.join([letter]*depth)

    def _str_dashgoid(self, dct, depth, single_or_double):
        """Return a string containing dashes (optional) and GO ID."""
        return "{DASHES} {ID}".format(
            DASHES=self._str_dash(depth, single_or_double) if self.indent else "",
            ID=self.nm2prtfmt['ID'].format(**dct))


# Copyright (C) 2016-2018, DV Klopfenstein, H Tang. All rights reserved.
