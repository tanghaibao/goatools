"""Get descendant/parent counts for all GO terms in a GODag and broad L0 and L1 terms."""

from __future__ import print_function

__copyright__ = "Copyright (C) 2016-present, DV Klopfenstein, H Tang, All rights reserved."
__author__ = "DV Klopfenstein"

import collections as cx
from itertools import chain
from goatools.godag.go_tasks import get_go2ancestors
from goatools.godag.go_tasks import get_go2descendants
from goatools.gosubdag.go_tasks import get_goobjs_altgo2goobj
from goatools.gosubdag.go_tasks import add_alt_goids


class CountRelativesInit:
    """Get descendant/parent counts for all GO terms in a GODag and broad L0 and L1 terms."""


    def __init__(self, go2obj, relationships, dcnt, go2letter):
        # Subset go2obj contains only items needed by go_sources
        self.go2obj = go2obj
        self.relationships = relationships
        self.dcnt = dcnt
        self.go2letter = go2letter
        # Ex: set(['part_of', 'regulates', 'negatively_regulates', 'positively_regulates'])
        _goobjs, _altgo2goobj = get_goobjs_altgo2goobj(self.go2obj)
        _r0 = not relationships  # True if not using relationships
        self.go2descendants = get_go2descendants(_goobjs, relationships)
        self.go2ancestors = get_go2ancestors(_goobjs, relationships)
        self.go2dcnt = cx.Counter({go: len(p) for go, p in self.go2descendants.items()})
        add_alt_goids(self.go2ancestors, _altgo2goobj)
        add_alt_goids(self.go2descendants, _altgo2goobj)
        add_alt_goids(self.go2dcnt, _altgo2goobj)
        # print('INIT CountRelativesInit', self.relationships)

    ## def get_relationship_dicts(self):
    ##     """Given GO DAG relationships, return summaries per GO ID."""
    ##     if not self.relationships:
    ##         return None
    ##     for goid, goobj in self.go2obj.items():
    ##         for reltyp, relset in goobj.relationship.items():
    ##             relfwd_goids = set(o.id for o in relset)
    ##             # for relfwd_goid in relfwd_goids:
    ##             #     assert relfwd_goid in self.go2obj, "{GO} {REL} NOT FOUND {GO_R}".format(
    ##             #         GO=goid, REL=reltyp, GO_R=relfwd_goid)
    ##             # print("CountRelativesInit RELLLLS", goid, goobj.id, reltyp, relfwd_goids)

    def get_goone2ntletter(self, go2dcnt, depth2goobjs):
        """Assign letters to depth-01 GO terms ordered using descendants cnt."""
        # 1. Group level-01/depth-01 GO terms by namespace
        ns2dcntgoobj = cx.defaultdict(list)
        for goobj in depth2goobjs[1]:
            dcnt = go2dcnt[goobj.id]
            ns2dcntgoobj[goobj.namespace].append((dcnt, goobj))
        # 2. Assign letters to level-01/depth-01 GO terms
        go2nt = {}
        ntobj = cx.namedtuple("NtGoLetters", "D1 dcnt goobj")
        _go2abc = self.go2letter
        letters = list(chain(range(ord('A'), ord('Z') + 1), range(ord('a'), ord('z') + 1)))
        for list_dcnt_goobj in ns2dcntgoobj.values():
            letter_idx = 0
            for dcnt, goobj in sorted(list_dcnt_goobj, key=lambda t: t[0], reverse=True):
                letter = chr(letters[letter_idx]) if _go2abc is None else _go2abc.get(goobj.id, '')
                go2nt[goobj.id] = ntobj._make([letter, dcnt, goobj])
                letter_idx += 1
        return go2nt

    @staticmethod
    def get_depth2goobjs(go2obj, max_depth=2):
        """Init depth2goobjs using list sorted by depth, get level-00/01 GO terms."""
        depth2goobjs = {d:list() for d in range(max_depth+1)}
        goid_seen = set()
        for _, goobj in sorted(go2obj.items(), key=lambda t: t[1].depth):
            # Save depth-00, depth-01, depth-02
            if goobj.depth > max_depth:
                break
            goid = goobj.id
            if not goobj.is_obsolete and goid not in goid_seen:
                depth2goobjs[goobj.depth].append(goobj)
                goid_seen.add(goid)
        return depth2goobjs


# Copyright (C) 2016-present, DV Klopfenstein, H Tang, All rights reserved.
