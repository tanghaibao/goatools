"""Create Typedef object."""

__copyright__ = "Copyright (C) 2015-2018, DV Klopfenstein, H Tang, All rights reserved."
__author__ = "DV Klopfenstein"


# pylint: disable=too-few-public-methods
class TypeDef(object):
    """
        TypeDef term.
    """

    def __init__(self):
        self.id = ""                # GO:NNNNNNN **DEPRECATED** reserved in Python
        self.item_id = ""           # GO:NNNNNNN (replaces deprecated "id")
        self.name = ""              # description
        self.namespace = ""         # external
        self.transitive_over = []   # List of other typedefs
        self.inverse_of = ""        # Name of inverse typedef.

    def __str__(self):
        ret = []
        ret.append("Typedef - {} ({}):".format(self.item_id, self.name))
        ret.append("  Inverse of: {}".format(self.inverse_of
                                             if self.inverse_of else "None"))
        if self.transitive_over:
            ret.append("  Transitive over:")
            for txo in self.transitive_over:
                ret.append("    - {}".format(txo))
        return "\n".join(ret)


def add_to_typedef(typedef_curr, obo_line):
    """Add new fields to the current typedef."""
    if obo_line[:4] == "id: ":
        assert not typedef_curr.item_id
        item_id = obo_line[4:]
        typedef_curr.item_id = item_id
    elif obo_line[:6] == "name: ":
        assert not typedef_curr.name
        typedef_curr.name = obo_line[6:]
    elif obo_line[:11] == "namespace: ":
        assert not typedef_curr.namespace
        typedef_curr.namespace = obo_line[11:]
    elif obo_line[17:] == "transitive_over: ":
        field_value = obo_line[17:].split('!')[0].rstrip()
        typedef_curr.transitive_over.append(field_value)
    elif obo_line[12:] == "inverse_of":
        assert not typedef_curr.inverse_of
        field_value = obo_line[12:].split('!')[0].rstrip()
        typedef_curr.inverse_of = field_value


# Copyright (C) 2015-2018, DV Klopfenstein, H Tang, All rights reserved.
