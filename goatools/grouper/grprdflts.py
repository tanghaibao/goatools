"""Given user GO ids and parent terms, group user GO ids under one parent term.

   Given a group of GO ids with one or more higher-level grouping terms, group
   each user GO id under the most descriptive parent GO term.

   Each GO id may have more than one parent.  One of the parent(s) is chosen
   to best represent the user GO id's function. The choice of parent is made by
   regarding how close the parent GO id is to the bottom of its hierarchy.

   The estimation of how close a GO term is to "the bottom" of its GO hierarchy
   is estimated using the number of total Go term descendent counts below
   that term.
"""

__copyright__ = "Copyright (C) 2016-2018, DV Klopfenstein, H Tang, All rights reserved."
__author__ = "DV Klopfenstein"

import sys
from goatools.base import get_godag
from goatools.gosubdag.gosubdag import GoSubDag


class GrouperDflts(object):
    """Holds objects that we would like to load and initialize once.

       Files used for grouping GO IDs:

           http://geneontology.org/ontology/go-basic.obo
           http://geneontology.org/ontology/subsets/goslim_generic.obo

    """

    def __init__(self, gosubdag=None, goslim_filename="goslim_generic.obo", hdrgos=None):
        self.gosubdag = self.get_gosubdag(gosubdag)
        _dagslim = get_godag(goslim_filename, prt=sys.stdout, loading_bar=False)
        self.ver_goslims = _dagslim.version
        self.goslims = self._init_goslims(_dagslim)
        self.hdrgos_dflt = self._init_hdrgos() if hdrgos is None else hdrgos  # goid set

    def _init_hdrgos(self):
        """Return GO IDs used as the default for the high grouping GO IDs."""
        # Get all GO terms that are at depth-00 or depth-01
        hdrgos = self.get_gos_d0d1()
        hdrgos |= self.goslims
        # self.gosubdag.prt_goids(hdrgos)
        return hdrgos

    def _init_goslims(self, dagslim):
        """Get GO IDs in GO slims."""
        go2obj_main = self.gosubdag.go2obj
        go2obj_slim = {go for go, o in dagslim.items() if go in go2obj_main}
        if self.gosubdag.relationships:
            return self._get_goslimids_norel(go2obj_slim)
        return set(dagslim.keys())

    def get_gos_d0d1(self):
        """Return GO IDs whose depth is 0 (BP, MF, CC) or depth is 1."""
        return set([o.id for d in [0, 1] for o in self.gosubdag.rcntobj.depth2goobjs.get(d)])

    def _get_goslimids_norel(self, dagslim):
        """Get all GO slim GO IDs that do not have a relationship."""
        go_slims = set()
        go2obj = self.gosubdag.go2obj
        for goid in dagslim:
            goobj = go2obj[goid]
            if not goobj.relationship:
                go_slims.add(goobj.id)
        return go_slims

    @staticmethod
    def get_gosubdag(gosubdag=None):
        """Gets a GoSubDag initialized for use by a Grouper object."""
        if gosubdag is not None:
            if gosubdag.rcntobj is not None:
                return gosubdag
            else:
                gosubdag.init_auxobjs()
                return gosubdag
        else:
            go2obj = get_godag()
            return GoSubDag(None, go2obj, rcntobj=True)


# Copyright (C) 2016-2018, DV Klopfenstein, H Tang, All rights reserved.
