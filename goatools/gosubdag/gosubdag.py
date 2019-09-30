"""Manages a user-specified subset of a GO DAG."""

from __future__ import print_function

__copyright__ = "Copyright (C) 2016-2019, DV Klopfenstein, H Tang, All rights reserved."
__author__ = "DV Klopfenstein"

import sys
import collections as cx
# import timeit
# from goatools.godag.prttime import prt_hms
from goatools.gosubdag.gosubdag_init import InitGOs
from goatools.gosubdag.gosubdag_init import InitFields
from goatools.gosubdag.go_tasks import chk_goids


class GoSubDag(object):
    """Manages a user-specified subset of a GO DAG."""

    def __init__(self, go_sources, go2obj, relationships=None, **kws):
        # kws _Init: rcntobj
        # tic = timeit.default_timer()
        _ini = InitGOs(go_sources, go2obj, relationships, **kws)
        self.go_sources = _ini.go_sources # set(go_sources)
        self.go2obj = _ini.go2obj # go2obj # Initialized with goobjs corresponding to go_sources
        self.relationships = _ini.relationships
        # tic = prt_hms(tic, "GoSubDag: InitGOs")
        # GO IDs to total count of all descendants: Init to None or CountRelatives object
        _fld = InitFields(_ini, **kws)
        self.rcntobj = _fld.get_rcntobj()  # None or CountRelatives object
        self.prt_attr = {
            'flds':_fld.prt_flds,           # namedtuple fields in go2nt
            'fmt':_fld.get_prt_fmt(False),  # GO:NNNNNNN   No indication if an alternate GO ID
            'fmta':_fld.get_prt_fmt(True)}  # GO:NNNNNNNa  'a' indicates if an alternate GO ID
        ### tic = _rpt_hms(tic, "GoSubDag: Create GoDepth1Letters")
        self.go2nt = _fld.get_go2nt_all(self.rcntobj)
        ### tic = _rpt_hms(tic0, "GoSubDag: total")
        prt = kws.get('prt', sys.stdout)
        if prt:
            self.prt_objdesc(prt)

    def prt_goids(self, goids=None, prtfmt=None, sortby=True, prt=sys.stdout):
        """Given GO IDs, print decriptive info about each GO Term."""
        if goids is None:
            goids = self.go_sources
        nts = self.get_nts(goids, sortby)
        if prtfmt is None:
            prtfmt = self.prt_attr['fmta']
        for ntgo in nts:
            key2val = ntgo._asdict()
            prt.write("{GO}\n".format(GO=prtfmt.format(**key2val)))
        return nts

    def get_nts(self, goids=None, sortby=None):
        """Given GO IDs, get a list of namedtuples."""
        nts = []
        # User GO IDs
        if goids is None:
            goids = self.go_sources
        else:
            chk_goids(goids, "GoSubDag::get_nts")
        if goids:
            ntobj = cx.namedtuple("NtGo", " ".join(self.prt_attr['flds']))
            go2nt = self.get_go2nt(goids)
            for goid, ntgo in self._get_sorted_go2nt(go2nt, sortby):
                assert ntgo is not None, "{GO} NOT IN go2nt".format(GO=goid)
                if goid == ntgo.GO:
                    nts.append(ntgo)
                else:
                    fld2vals = ntgo._asdict()
                    fld2vals['GO'] = goid
                    nts.append(ntobj(**fld2vals))
        return nts

    def _get_sorted_go2nt(self, go2nt, sortby):
        """Return sorted list of tuples."""
        if sortby is True:
            _fnc = self.get_fncsortnt()
            return sorted(go2nt.items(), key=lambda t: _fnc(t[1]))
        if sortby:
            return sorted(go2nt.items(), key=lambda t: sortby(t[1]))
        return go2nt.items()

    def get_fncsortnt(self):
        """Return sorted list of tuples."""
        if 'dcnt' in self.prt_attr['flds']:
            if 'D1' in self.prt_attr['flds']:
                return lambda ntgo: [ntgo.NS, ntgo.depth, -1*ntgo.dcnt, ntgo.D1, ntgo.alt]
            else:
                return lambda ntgo: [ntgo.NS, ntgo.depth, -1*ntgo.dcnt, ntgo.alt]
        else:
            return lambda ntgo: [ntgo.NS, -1*ntgo.depth, ntgo.alt]

    def get_go2nt(self, goids):
        """Return dict of GO ID as key and GO object information in namedtuple."""
        get_nt = self.go2nt
        goids_present = set(goids).intersection(self.go2obj)
        if len(goids_present) != len(goids):
            print("GO IDs NOT FOUND IN DAG: {GOs}".format(
                GOs=" ".join(set(goids).difference(goids_present))))
        return {g:get_nt[g] for g in goids_present}

    def get_go2obj(self, goids):
        """Return a go2obj dict for just the user goids."""
        go2obj = self.go2obj
        return {go:go2obj[go] for go in goids}

    def get_vals(self, field, goids=None):
        """Return a go2obj dict for just the user goids."""
        go2nt = self.go2nt
        if goids is None:
            goids = set(go2nt)
        return [getattr(go2nt[go], field) for go in goids]

    def get_key_goids(self, goids):
        """Given GO IDs, return key GO IDs."""
        go2obj = self.go2obj
        return set(go2obj[go].id for go in goids)

    def get_ns2goids(self, goids):
        """Group GO IDs by namespace."""
        ns2goids = cx.defaultdict(set)
        go2nt = self.go2nt
        for goid in goids:
            ns2goids[go2nt[goid].NS].add(goid)
        return {ns:gos for ns, gos in ns2goids.items()}

    def prt_objdesc(self, prt):
        """Return description of this GoSubDag object."""
        txt = "INITIALIZING GoSubDag: {N:3} sources in {M:3} GOs rcnt({R}). {A} alt GO IDs\n"
        alt2obj = {go:o for go, o in self.go2obj.items() if go != o.id}
        prt.write(txt.format(
            N=len(self.go_sources),
            M=len(self.go2obj),
            R=self.rcntobj is not None,
            A=len(alt2obj)))
        prt.write("             GoSubDag: namedtuple fields: {FLDS}\n".format(
            FLDS=" ".join(self.prt_attr['flds'])))
        prt.write("             GoSubDag: relationships: {RELS}\n".format(RELS=self.relationships))


# Copyright (C) 2016-2019, DV Klopfenstein, H Tang, All rights reserved.
