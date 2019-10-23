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

import sys
import collections as cx
from goatools.godag.consts import NAMESPACE2NS
from goatools.grouper.grprobj_init import GrouperInit

__copyright__ = "Copyright (C) 2016-2019, DV Klopfenstein, H Tang, All rights reserved."
__author__ = "DV Klopfenstein"


class Grouper(object):
    """Groups the user GO ids under other GO IDs acting as headers for the GO groups."""

    fmtsum = ("{GO_DESC} GOs({GOs:6,} in {SECs:2} sections, "
              "{UNGRP:>3} {undesc}) {ACTION} {FILE}\n")

    def __init__(self, grpname, goids, hdrobj, gosubdag, **kws):
        # print("INITIALIZING Grouper")
        # Data members read
        self.grpname = grpname
        self.hdrobj = hdrobj  # Contains all possible hdrgos, not just ones used
        self.gosubdag = gosubdag
        assert self.gosubdag.rcntobj is not None
        # _ini = GrouperInit(grpname, goids, hdrobj, gosubdag, kws.get('fnc_most_specific', 'dcnt'))
        _ini = GrouperInit(goids, self, kws.get('fnc_most_specific', 'dcnt'))
        self.usrgos = _ini.usrgos
        # Initialize: hdrgo2usrgos hdrgo_is_usrgo
        #   * hdrgo2usrgos: User GO IDs, grouped under high GO IDs (grouped, but not sorted)
        self.hdrgo2usrgos = _ini.hdrgo2usrgos
        self.hdrgo_is_usrgo = _ini.hdrgo_is_usrgo  # set of GO IDs -> both headers/user GO IDs
        # User GO IDs and their corresponding high GO IDs (not grouped or sorted)
        self.go2nt = _ini.get_go2nt(kws.get('go2nt', None))

    def get_usrgos_w_parents(self, hdrgos, usrgos_all=None):
        """Get usrgos w/parents in hdrgos, even if usrgos did not get grouped under hdrgos."""
        usrgos = set()
        _go2parents = self.gosubdag.rcntobj.go2ancestors
        if usrgos_all is None:
            usrgos_all = self.usrgos
        for usrgo in usrgos_all:
            all_usrgo_parents = _go2parents.get(usrgo)
            sel_usrgo_parents = all_usrgo_parents.intersection(hdrgos)
            if sel_usrgo_parents:
                usrgos.add(usrgo)
        return usrgos

    def get_sections_2d(self):
        """Get 2-D list of sections and hdrgos sets actually used in grouping."""
        sections_hdrgos_act = []
        hdrgos_act_all = self.get_hdrgos()  # Header GOs actually used to group
        hdrgos_act_secs = set()
        if self.hdrobj.sections:
            for section_name, hdrgos_all_lst in self.hdrobj.sections:
                # print("GGGGGGGGGGGGGGGGG {N:3} {NAME}".format(N=len(hdrgos_all_lst), NAME=section_name))
                hdrgos_all_set = set(hdrgos_all_lst)
                hdrgos_act_set = hdrgos_all_set.intersection(hdrgos_act_all)
                if hdrgos_act_set:
                    hdrgos_act_secs |= hdrgos_act_set
                    # Use original order of header GOs found in sections
                    hdrgos_act_lst = []
                    hdrgos_act_ctr = cx.Counter()
                    for hdrgo_p in hdrgos_all_lst: # Header GO that may or may not be used.
                        if hdrgo_p in hdrgos_act_set and hdrgos_act_ctr[hdrgo_p] == 0:
                            hdrgos_act_lst.append(hdrgo_p)
                        hdrgos_act_ctr[hdrgo_p] += 1
                    sections_hdrgos_act.append((section_name, hdrgos_act_lst))
            # print(">>>>>>>>>>>>>>> hdrgos_act_all {N:3}".format(N=len(hdrgos_act_all)))
            # print(">>>>>>>>>>>>>>> hdrgos_act_secs {N:3}".format(N=len(hdrgos_act_secs)))
            hdrgos_act_rem = hdrgos_act_all.difference(hdrgos_act_secs)
            if hdrgos_act_rem:
                # print("RRRRRRRRRRR {N:3}".format(N=len(hdrgos_act_rem)))
                sections_hdrgos_act.append((self.hdrobj.secdflt, hdrgos_act_rem))
        else:
            sections_hdrgos_act.append((self.hdrobj.secdflt, hdrgos_act_all))
        return sections_hdrgos_act

    def get_usrgos_g_section(self, section=None):
        """Get usrgos in a requested section."""
        if section is None:
            section = self.hdrobj.secdflt
        if section is True:
            return self.usrgos
        # Get dict of sections and hdrgos actually used in grouping
        section2hdrgos = cx.OrderedDict(self.get_sections_2d())
        hdrgos_lst = section2hdrgos.get(section, None)
        if hdrgos_lst is not None:
            hdrgos_set = set(hdrgos_lst)
            hdrgos_u = hdrgos_set.intersection(self.hdrgo_is_usrgo)
            hdrgos_h = hdrgos_set.intersection(self.hdrgo2usrgos.keys())
            usrgos = set([u for h in hdrgos_h for u in self.hdrgo2usrgos.get(h)])
            usrgos |= hdrgos_u
            return usrgos
        return set()

    def get_section2usrnts(self):
        """Get dict section2usrnts."""
        sec_nts = []
        for section_name, _ in self.get_sections_2d():
            usrgos = self.get_usrgos_g_section(section_name)
            sec_nts.append((section_name, [self.go2nt.get(u) for u in usrgos]))
        return cx.OrderedDict(sec_nts)

    def get_section2items(self, itemkey):
        """Collect all items into a single set per section."""
        sec_items = []
        section2usrnts = self.get_section2usrnts()
        for section, usrnts in section2usrnts.items():
            items = set([e for nt in usrnts for e in getattr(nt, itemkey, set())])
            sec_items.append((section, items))
        return cx.OrderedDict(sec_items)

    def get_hdrgos_g_usrgos(self, usrgos):
        """Return hdrgos which contain the usrgos."""
        hdrgos_for_usrgos = set()
        hdrgos_all = self.get_hdrgos()
        usrgo2hdrgo = self.get_usrgo2hdrgo()
        for usrgo in usrgos:
            if usrgo in hdrgos_all:
                hdrgos_for_usrgos.add(usrgo)
                continue
            hdrgo_cur = usrgo2hdrgo.get(usrgo, None)
            if hdrgo_cur is not None:
                hdrgos_for_usrgos.add(hdrgo_cur)
        return hdrgos_for_usrgos

    def get_section_hdrgos_nts(self, sortby=None):
        """Get a flat list of sections and hdrgos actually used in grouping."""
        nts_all = []
        section_hdrgos_actual = self.get_sections_2d()
        flds_all = ['Section'] + self.gosubdag.prt_attr['flds']
        ntobj = cx.namedtuple("NtGoSec", " ".join(flds_all))
        flds_go = None
        if sortby is None:
            sortby = lambda nt: -1*nt.dcnt
        for section_name, hdrgos_actual in section_hdrgos_actual:
            nts_sec = []
            for hdrgo_nt in self.gosubdag.get_go2nt(hdrgos_actual).values():
                if flds_go is None:
                    flds_go = hdrgo_nt._fields
                key2val = {key:val for key, val in zip(flds_go, list(hdrgo_nt))}
                key2val['Section'] = section_name
                nts_sec.append(ntobj(**key2val))
            nts_all.extend(sorted(nts_sec, key=sortby))
        return nts_all

    def get_sections_2d_nts(self, sortby=None):
        """Get high GO IDs that are actually used to group current set of GO IDs."""
        sections_2d_nts = []
        for section_name, hdrgos_actual in self.get_sections_2d():
            hdrgo_nts = self.gosubdag.get_nts(hdrgos_actual, sortby=sortby)
            sections_2d_nts.append((section_name, hdrgo_nts))
        return sections_2d_nts

    def get_hdrgos(self):
        """Return high GO IDs that are actually used to group current set of GO IDs."""
        return set(self.hdrgo2usrgos.keys()).union(self.hdrgo_is_usrgo)

    def get_usrgos_g_hdrgos(self, hdrgos):
        """Return usrgos under provided hdrgos."""
        usrgos_all = set()
        if isinstance(hdrgos, str):
            hdrgos = [hdrgos]
        for hdrgo in hdrgos:
            usrgos_cur = self.hdrgo2usrgos.get(hdrgo, None)
            if usrgos_cur is not None:
                usrgos_all |= usrgos_cur
            if hdrgo in self.hdrgo_is_usrgo:
                usrgos_all.add(hdrgo)
        return usrgos_all

    def get_hdrgos_unplaced(self):
        """Get hdrgos which are not headers in sections."""
        return self.get_hdrgos().difference(self.hdrobj.get_section_hdrgos())

    def get_hdrgos_u0(self):
        """Return header GO IDs which ARE NOT user GO IDs."""
        return set(self.hdrgo2usrgos.keys()).difference(self.usrgos)

    def get_hdrgos_u1(self):
        """Return header GO IDs which ARE user GO IDs."""
        return self.hdrgo_is_usrgo

    def get_hdrgo2usrgos(self, hdrgos):
        """Return a subset of hdrgo2usrgos."""
        get_usrgos = self.hdrgo2usrgos.get
        hdrgos_actual = self.get_hdrgos().intersection(hdrgos)
        return {h:get_usrgos(h) for h in hdrgos_actual}

    def get_usrgo2hdrgo(self):
        """Return a dict with all user GO IDs as keys and their respective header GOs as values."""
        usrgo2hdrgo = {}
        for hdrgo, usrgos in self.hdrgo2usrgos.items():
            for usrgo in usrgos:
                assert usrgo not in usrgo2hdrgo
                usrgo2hdrgo[usrgo] = hdrgo
        # Add usrgos which are also a hdrgo and the GO group contains no other GO IDs
        for goid in self.hdrgo_is_usrgo:
            usrgo2hdrgo[goid] = goid
        assert len(self.usrgos) <= len(usrgo2hdrgo), \
            "USRGOS({U}) != USRGO2HDRGO({H}): {GOs}".format(
                U=len(self.usrgos),
                H=len(usrgo2hdrgo),
                GOs=self.usrgos.symmetric_difference(set(usrgo2hdrgo.keys())))
        return usrgo2hdrgo

    def get_go2sectiontxt(self):
        """Return a dict with actual header and user GO IDs as keys and their sections as values."""
        go2txt = {}
        _get_secs = self.hdrobj.get_sections
        hdrgo2sectxt = {h:" ".join(_get_secs(h)) for h in self.get_hdrgos()}
        usrgo2hdrgo = self.get_usrgo2hdrgo()
        for goid, ntgo in self.go2nt.items():
            hdrgo = ntgo.GO if ntgo.is_hdrgo else usrgo2hdrgo[ntgo.GO]
            go2txt[goid] = hdrgo2sectxt[hdrgo]
        return go2txt

    def get_usrgo2sections(self):
        """Return a dict with all user GO IDs as keys and their sections as values."""
        usrgo2sections = cx.defaultdict(set)
        usrgo2hdrgo = self.get_usrgo2hdrgo()
        get_sections = self.hdrobj.get_sections
        for usrgo, hdrgo in usrgo2hdrgo.items():
            sections = set(get_sections(hdrgo))
            usrgo2sections[usrgo] |= sections
        assert len(usrgo2sections) >= len(self.usrgos), \
            "uGOS({U}) != uGO2sections({H}): {GOs}".format(
                U=len(self.usrgos),
                H=len(usrgo2sections),
                GOs=self.usrgos.symmetric_difference(set(usrgo2sections.keys())))
        return usrgo2sections

    def get_fout_base(self, goid, name=None, pre="gogrp"):
        """Get filename for a group of GO IDs under a single header GO ID."""
        goobj = self.gosubdag.go2obj[goid]
        if name is None:
            name = self.grpname.replace(" ", "_")
        sections = "_".join(self.hdrobj.get_sections(goid))
        return "{PRE}_{BP}_{NAME}_{SEC}_{DSTR}_{D1s}_{GO}".format(
            PRE=pre,
            BP=NAMESPACE2NS[goobj.namespace],
            NAME=self._str_replace(name),
            SEC=self._str_replace(self._str_replace(sections)),
            GO=goid.replace(":", ""),
            DSTR=self._get_depthsr(goobj),
            D1s=self.gosubdag.go2nt[goobj.id].D1)

    def _get_depthsr(self, goobj):
        """Return DNN or RNN depending on if relationships are loaded."""
        if 'reldepth' in self.gosubdag.prt_attr['flds']:
            return "R{R:02}".format(R=goobj.reldepth)
        return "D{D:02}".format(D=goobj.depth)

    @staticmethod
    def _str_replace(txt):
        """Makes a small text amenable to being used in a filename."""
        txt = txt.replace(",", "")
        txt = txt.replace(" ", "_")
        txt = txt.replace(":", "")
        txt = txt.replace(".", "")
        txt = txt.replace("/", "")
        txt = txt.replace("", "")
        return txt

    def prt_summary(self, prt=sys.stdout):
        """Print summary of grouping/sorting run."""
        # Grouping summary
        fmtstr = "Grouped: {U:3,} User GOs, using {h:2,} of {H:,} Grouping GOs, for run: {NAME}\n"
        prt.write(fmtstr.format(
            NAME=self.grpname,
            U=len(self.usrgos),
            h=len(self.hdrobj.hdrgos.intersection(self.hdrgo2usrgos.keys())),
            H=self.hdrobj.num_hdrgos()))

# Copyright (C) 2016-2019, DV Klopfenstein, H Tang, All rights reserved.
