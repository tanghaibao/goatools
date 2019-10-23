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

from __future__ import print_function

import collections as cx
from goatools.nt_utils import get_dict_w_id2nts
from goatools.gosubdag.go_most_specific import get_most_specific_dcnt
from goatools.gosubdag.go_most_specific import get_most_specific_tinfo
from goatools.gosubdag.go_most_specific import get_most_specific_tinfo_dcnt
from goatools.grouper.utils import get_hdridx_flds

__copyright__ = "Copyright (C) 2016-2018, DV Klopfenstein, H Tang, All rights reserved."
__author__ = "DV Klopfenstein"


class GrouperInit:
    """Initialize Grouper object."""

    most_specific_fncs = {
        'dcnt': get_most_specific_dcnt,
        'tinfo': get_most_specific_tinfo,
        'tinfo_dcnt': get_most_specific_tinfo_dcnt}

    def __init__(self, goids, objpre, fnc_most_specific='dcnt'):
        # Data members read
        self.grpname = objpre.grpname
        self.gosubdag = objpre.gosubdag
        self.usrgos = self._init_usrgos(goids)
        self.hdrobj = objpre.hdrobj # Contains all possible hdrgos, not just ones used
        assert self.gosubdag.rcntobj is not None
        # Initialize: hdrgo2usrgos hdrgo_is_usrgo
        #   * hdrgo2usrgos: User GO IDs, grouped under high GO IDs (grouped, but not sorted)
        self.hdrgo2usrgos = None
        self.hdrgo_is_usrgo = None  # Will contain both main GO IDs and user-specified alt GO IDs
        self._init_h2us(fnc_most_specific)

    def _init_usrgos(self, goids):
        """Return user GO IDs which have GO Terms."""
        usrgos = set()
        goids_missing = set()
        _go2obj = self.gosubdag.go2obj
        for goid in goids:
            if goid in _go2obj:
                usrgos.add(goid)
            else:
                goids_missing.add(goid)
        if goids_missing:
            print("MISSING GO IDs: {GOs}".format(GOs=goids_missing))
            print("{N} of {M} GO IDs ARE MISSING".format(N=len(goids_missing), M=len(goids)))
        return usrgos

    def get_gos_all(self):
        """Return a flat list of all GO IDs in grouping object.

           All GO IDs:
               * header GO IDs that are not user GO IDs
               * user GO IDs that are under header GOs
               * user GO IDs that are header GOs in groups containing no other user GO IDs
        """
        gos_all = set()
        # Get:
        #   * Header GO IDs that are not user GO IDs
        #   * User GO IDs that are under header GOs
        for hdrgo, usrgos in self.hdrgo2usrgos.items():
            gos_all.add(hdrgo)
            gos_all |= usrgos
        # User GO IDs that are header GOs in groups containing no other user GO IDs
        gos_all |= self.hdrgo_is_usrgo
        assert gos_all == self.usrgos.union(set(self.hdrgo2usrgos.keys()))
        assert not self.usrgos.difference(gos_all), \
            "GROUPER ERROR: {GOs}".format(GOs=self.usrgos.difference(gos_all))
        return gos_all

    def _init_h2us(self, fnc_most_specific):
        """Given a set of user GO ids, return GO ids grouped under the "GO high" terms.

        Example of a grouped go list:

            gos = ['GO:0044464':[    # grp_term: D1 cell part
                       'GO:0005737', # child:    D3 cytoplasm
                       'GO:0048471', # child:    D4 perinuclear region of cytoplasm
                   'GO:0016020':[    # grp_term: D1 membrane
                       'GO:0098589', # child:    D2 membrane region
                       'GO:0005886', # child:    D2 plasma membrane
                  ]
        """
        # Header GO IDs are main. User GO IDs are as specified by the user
        hdrgo2usrgos = cx.defaultdict(set)
        # Contains user GO IDs which are also header GO IDs, plus user main GO if needed
        hdrgo_is_usrgo = set()
        _go2nt = self.gosubdag.go2nt
        objhi = GrouperInit.GetGoidHigh(self.gosubdag, self.hdrobj.hdrgos,
                                        self.most_specific_fncs[fnc_most_specific])
        for goid_usr in self.usrgos:
            goid_main = _go2nt[goid_usr].id
            # Add current GO ID to parents_all in case curr GO ID is a high GO.
            goid_high = objhi.get_goid_high(goid_main)
            # Don't add user GO ID if it is also the GO header
            if goid_main != goid_high:
                hdrgo2usrgos[goid_high].add(goid_usr)
            elif goid_high not in hdrgo2usrgos:
                hdrgo2usrgos[goid_high] = set()
            if goid_main == goid_high:
                hdrgo_is_usrgo.add(goid_main)
                if goid_main != goid_usr:
                    hdrgo_is_usrgo.add(goid_usr)
        # Initialize data members
        self.hdrgo2usrgos = hdrgo2usrgos
        self.hdrgo_is_usrgo = hdrgo_is_usrgo

    # pylint: disable=too-few-public-methods
    class GetGoidHigh:
        """Given a user GO ID, return the 'closest' header GO."""

        def __init__(self, gosubdag, gos_high, get_most_specific):
            self.go2parents = gosubdag.rcntobj.go2ancestors
            self.go2nt = gosubdag.go2nt
            self.gos_high = gos_high
            self.get_most_specific = get_most_specific

        def get_goid_high(self, goid_main):
            """Return the 'closest' GO header to the GO ID arg."""
            parents_all = {goid_main}
            if goid_main in self.go2parents:
                parents_all.update(self.go2parents[goid_main])
            parents_high = parents_all.intersection(self.gos_high)
            assert parents_high, "NO PARENTS {P} {H} {NT}".format(
                P=len(parents_all), H=len(self.gos_high), NT=goid_main)
            return self.get_most_specific(parents_high, self.go2nt)


    # --- Initialize go2nt. Namedtuple fields may be used in sortby lambda functions
    def get_go2nt(self, usr_go2nt):
        """Combine user namedtuple fields, GO object fields, and format_txt."""
        gos_all = self.get_gos_all()
        # Minimum set of namedtuple fields available for use with Sorter on grouped GO IDs
        prt_flds_all = get_hdridx_flds() + self.gosubdag.prt_attr['flds']
        if not usr_go2nt:
            return self.__init_go2nt_dflt(gos_all, prt_flds_all)
        usr_nt_flds = next(iter(usr_go2nt.values()))._fields
        # If user namedtuple already contains all fields available, then return usr_go2nt
        if not set(prt_flds_all).difference(usr_nt_flds):
            return self._init_go2nt_aug(usr_go2nt)
        # Otherwise, combine user fields and default Sorter fields
        return self.__init_go2nt_w_usr(gos_all, usr_go2nt, prt_flds_all)

    def __init_go2nt_dflt(self, gos_all, prt_flds_all):
        """Combine GO object fields and format_txt."""
        go2nts = [self.gosubdag.go2nt, self._get_go2nthdridx(gos_all)]
        go2nt = get_dict_w_id2nts(gos_all, go2nts, prt_flds_all)
        return self._init_go2nt_aug(go2nt)

    def __init_go2nt_w_usr(self, gos_all, usr_go2nt, prt_flds_all):
        """Combine GO object fields and format_txt."""
        assert usr_go2nt, "go2nt HAS NO ELEMENTS"
        from goatools.nt_utils import get_unique_fields
        go2nts = [usr_go2nt, self.gosubdag.go2nt, self._get_go2nthdridx(gos_all)]
        usr_nt_flds = next(iter(usr_go2nt.values()))._fields # Get any single value from a dict
        flds = get_unique_fields([usr_nt_flds, prt_flds_all])
        go2nt = get_dict_w_id2nts(gos_all, go2nts, flds)
        return self._init_go2nt_aug(go2nt)

    def _init_go2nt_aug(self, go2nt):
        """Augment go2nt with GO ID key to account for alt GO IDs."""
        go2obj = self.gosubdag.go2obj
        # Get alt GO IDs
        go2nt_aug = {}
        # NOW
        for goid_usr, nt_usr in go2nt.items():
            goobj = go2obj[goid_usr]
            if goobj.alt_ids:
                alts = set(goobj.alt_ids)
                alts.add(goobj.id)
                for goid_alt in alts:
                    if goid_alt not in go2nt:
                        go2nt_aug[goid_alt] = nt_usr
        # WAS
        # Add alt GO IDs to go2nt
        for goid, gont in go2nt_aug.items():
            go2nt[goid] = gont
        return go2nt

    def _get_go2nthdridx(self, gos_all):
        """Get GO IDs header index for each user GO ID and corresponding parent GO IDs."""
        go2nthdridx = {}
        # NtHdrIdx Namedtuple fields:
        #   * format_txt: Used to determine the format when writing Excel cells
        #   * hdr_idx: Value printed in an Excel cell
        # shortcuts
        obj = GrouperInit.NtMaker(self)
        # Create go2nthdridx
        for goid in gos_all:
            go2nthdridx[goid] = obj.get_nt(goid)
        return go2nthdridx

    class NtMaker:
        """Make namedtuples for GO IDs in grouper."""

        ntobj = cx.namedtuple("NtHdrIdx", " ".join(get_hdridx_flds()))

        def __init__(self, obj):
            self.grpname = obj.grpname
            self.usrgos = obj.usrgos
            self.hdrgos = obj.hdrobj.hdrgos
            ## assert "GO:0008150" in self.hdrgos
            self.go2obj = obj.gosubdag.go2obj
            self.hdrgo2usrgos = obj.hdrgo2usrgos
            self.hdrgo_is_usrgo = obj.hdrgo_is_usrgo

        def get_nt(self, goid_user):
            """Get Grouper namedtuple for user GO ID."""
            goid_main = self.go2obj[goid_user].id
            goid_in_hdrgos = goid_main in self.hdrgo2usrgos
            goid_in_usrgos = goid_user in self.hdrgo_is_usrgo
            # format_txt = int(goid_in_hdrgos or goobj.id in self.hdrgos)
            format_txt = int(goid_in_hdrgos)
            # namedtuple grouping fields
            hdr1usr01 = self._get_hdr1usr01(goid_in_hdrgos, goid_in_usrgos)
            return self.ntobj(
                format_txt=format_txt,
                hdr_idx=format_txt,
                is_hdrgo=goid_in_hdrgos,
                is_usrgo=goid_in_usrgos,
                num_usrgos=self._get_num_usrgos(goid_user, goid_in_hdrgos, goid_in_usrgos),
                hdr1usr01=hdr1usr01)

        def _get_num_usrgos(self, goid_main, goid_in_hdrgos, goid_in_usrgos):
            """Get the number of user GO IDs under a header GO ID."""
            if not goid_in_hdrgos:
                return "."
            num_goids = len(self.hdrgo2usrgos[goid_main]) + int(goid_in_usrgos)
            assert num_goids != 0, "{NAME} MAIN({GO}) num_goids({N})\n{HDRUSR}".format(
                NAME=self.grpname, GO=goid_main, N=num_goids, HDRUSR=" ".join(sorted(self.hdrgos)))
            return num_goids

        @staticmethod
        def _get_hdr1usr01(goid_in_hdrgos, goid_in_usrgos):
            """Get string indicating if GO is also a header GO."""
            if goid_in_hdrgos:
                return "**" if goid_in_usrgos else "*"
            return ""


# Copyright (C) 2016-2018, DV Klopfenstein, H Tang, All rights reserved.
