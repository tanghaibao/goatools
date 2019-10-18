"""Get descendant/parent counts for all GO terms in a GODag and broad L0 and L1 terms."""

__copyright__ = "Copyright (C) 2016-present, DV Klopfenstein, H Tang, All rights reserved."
__author__ = "DV Klopfenstein"

from goatools.gosubdag.godag_rcnt_init import CountRelativesInit


class CountRelatives:
    """Get descendant/parent counts for all GO terms in a GODag and broad L0 and L1 terms."""

    def __init__(self, go2obj, relationships=None, dcnt=True, go2letter=None):
        # print("INITIALIZING CountRelatives")
        # Subset go2obj contains only items needed by go_sources
        self.go2obj = go2obj
        # Count of total number of descendants for each GO term
        _ini = CountRelativesInit(go2obj, relationships, dcnt, go2letter)
        self.go2descendants = _ini.go2descendants  # GO IDs
        # Used by: Semantic, Grouper
        self.go2parents = _ini.go2ancestors    # Will be DEPRECATED: renamed to go2ancestors
        self.go2ancestors = _ini.go2ancestors
        self.go2dcnt = _ini.go2dcnt
        # self.go_relationships = _ini.get_relationship_dicts()
        # Top (depth-00) terms (BP, MF, CC) and depth-01 terms
        self.depth2goobjs = _ini.get_depth2goobjs(go2obj)
        self.gos_depth1 = set(goobj.id for goobj in self.depth2goobjs[1])
        self.goone2ntletter = _ini.get_goone2ntletter(self.go2dcnt, self.depth2goobjs)

    def get_parents_letters(self, goobj):
        """Get the letters representing all parent terms which are depth-01 GO terms."""
        if goobj.id in self.go2ancestors:
            parents_all = set.union(self.go2ancestors[goobj.id])
            parents_all.add(goobj.id)
            # print "{}({}) D{:02}".format(goobj.id, goobj.name, goobj.depth), parents_all
            parents_d1 = parents_all.intersection(self.gos_depth1)
            return [self.goone2ntletter[g].D1 for g in parents_d1]
        return []

    def get_d1str(self, goobj, reverse=False):
        """Get D1-string representing all parent terms which are depth-01 GO terms."""
        return "".join(sorted(self.get_parents_letters(goobj), reverse=reverse))


# Copyright (C) 2016-present, DV Klopfenstein, H Tang, All rights reserved.
