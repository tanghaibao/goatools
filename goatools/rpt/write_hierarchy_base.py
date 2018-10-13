"""Print a GO term's lower-level hierarchy."""

__copyright__ = "Copyright (C) 2016-2018, DV Klopfenstein, H Tang. All rights reserved."
__author__ = "DV Klopfenstein"

import sys


# pylint: disable=too-many-instance-attributes,too-few-public-methods
class WrHierPrt(object):
    """Print the hierarchy of a set of objects w/attrs=[item_id, children]."""

    fmt_dashes = '{DASHES} {ID}'

    def __init__(self, id2obj, id2nt, cfg, prt=sys.stdout):
        self.id2obj = id2obj  # Contains children (and parents)
        self.id2nt = id2nt    # Contains fields for printing (optional)
        self.do_prtfmt = self.id2nt is not None
        self.nm2prtfmt = cfg['name2prtfmt'] if self.do_prtfmt else None
        self.max_indent = cfg['max_indent']
        self.include_only = cfg['include_only']
        self.item_marks = self._init_item_marks(cfg.get('item_marks'))
        self.mark_dflt = cfg.get('mark_default', ' ')
        self.concise_prt = cfg.get('concise_prt', False)
        self.indent = cfg.get('indent', True)
        self.space_branches = cfg.get('space_branches', False)
        self.sortby = cfg.get('sortby')
        # vars
        self.prt = prt
        self.items_printed = set()
        self.items_list = []
        self.dash_len = cfg['dash_len'] + cfg.get('id_len', 10) + 2

    def prt_hier_rec(self, item_id, depth=1):
        """Write hierarchy for a GO Term record and all GO IDs down to the leaf level."""
        # Shortens hierarchy report by only printing the hierarchy
        # for the sub-set of user-specified GO terms which are connected.
        if self.include_only and item_id not in self.include_only:
            return

        obj = self.id2obj[item_id]
        # Optionally space the branches for readability
        if self.space_branches:
            if depth == 1 and obj.children:
                self.prt.write("\n")
        # Print marks if provided
        if self.item_marks:
            self.prt.write('{MARK} '.format(
                MARK=self.item_marks.get(item_id, self.mark_dflt)))

        no_repeat = self.concise_prt and item_id in self.items_printed
        # Print content
        dashes = self._str_dash(depth, no_repeat, obj)
        if self.do_prtfmt:
            self._prtfmt(item_id, dashes)
        else:
            self._prtstr(obj, dashes)
        self.items_printed.add(item_id)
        self.items_list.append(item_id)
        # Do not print hierarchy below this turn if it has already been printed
        if no_repeat:
            return
        depth += 1
        if self.max_indent is not None and depth > self.max_indent:
            return
        children = obj.children if self.sortby is None else sorted(obj.children, key=self.sortby)
        for child in children:
            self.prt_hier_rec(child.item_id, depth)

    def _prtfmt(self, item_id, dashes):
        """Print object information using a namedtuple and a format pattern."""
        ntprt = self.id2nt[item_id]
        dct = ntprt._asdict()
        self.prt.write('{DASHES:{N}}'.format(
            DASHES=self.fmt_dashes.format(DASHES=dashes, ID=self.nm2prtfmt['ID'].format(**dct)),
            N=self.dash_len))
        self.prt.write("{INFO}\n".format(INFO=self.nm2prtfmt['ITEM'].format(**dct)))

    def _prtstr(self, obj, dashes):
        """Print object information using a namedtuple and a format pattern."""
        self.prt.write('{DASHES:{N}}'.format(
            DASHES=self.fmt_dashes.format(DASHES=dashes, ID=obj.item_id),
            N=self.dash_len))
        self.prt.write("{INFO}\n".format(INFO=str(obj)))

    def _str_dash(self, depth, no_repeat, obj):
        """Return a string containing dashes (optional) and GO ID."""
        if self.indent:
            # '-' is default character indicating hierarchy level
            # '=' is used to indicate a hierarchical path printed in detail previously.
            single_or_double = not no_repeat or not obj.children
            letter = '-' if single_or_double else '='
            return ''.join([letter]*depth)
        return ""

    @staticmethod
    def _init_item_marks(item_marks):
        """Initialize the makred item dict."""
        if isinstance(item_marks, dict):
            return item_marks
        if item_marks:
            return {item_id:'>' for item_id in item_marks}


# Copyright (C) 2016-2018, DV Klopfenstein, H Tang. All rights reserved.
