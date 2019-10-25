"""Manages a user-specified subset of a GO DAG."""

from __future__ import print_function

__copyright__ = "Copyright (C) 2016-2020, DV Klopfenstein, H Tang, All rights reserved."
__author__ = "DV Klopfenstein"

import sys
import collections as cx
import math
from goatools.godag.consts import NAMESPACE2NS
from goatools.godag.relationship_str import RelationshipStr
from goatools.godag.go_tasks import CurNHigher
from goatools.gosubdag.godag_rcnt import CountRelatives
from goatools.gosubdag.go_tasks import get_leaf_children
from goatools.gosubdag.utils import get_kwargs


# pylint: disable=too-few-public-methods
class InitGOs:
    """Initialize GoSubDab."""

    # Add additional GO IDs if used in user tasks
    kws_aux_gos = set(['go2color'])

    def __init__(self, go_sources, go2obj, relationships=False, **kws):
        # kws: go2color, children
        self.kws = kws
        # Process: rcntobj tcntobj go2nt relationships
        self.go2obj_orig = go2obj
        if relationships:
            assert hasattr(next(iter(go2obj.values())), 'relationship'), "NO DAG RELATIONSHIPS"
        # Init go2obj and go_sources
        self.go2obj = None
        self.go_sources = None
        self._init_gos(go_sources, relationships)
        # Using reduced go2obj, init relationships
        self.relationships = self._init_relationships(relationships)  # set of relationship types

    def _init_relationships(self, relationships_arg):
        """Return a set of relationships found in all subset GO Terms."""
        if relationships_arg:
            relationships_all = self._get_all_relationships()
            if relationships_arg is True:
                return relationships_all
            return relationships_all.intersection(relationships_arg)
        return set()

    def _get_all_relationships(self):
        """Return all relationships seen in GO Dag subset."""
        relationships_all = set()
        for goterm in self.go2obj.values():
            if goterm.relationship:
                relationships_all.update(goterm.relationship)
            if goterm.relationship_rev:
                relationships_all.update(goterm.relationship_rev)
        return relationships_all

    def _init_gos(self, go_sources_arg, relationships_arg):
        """Initialize GO sources."""
        # No GO sources provided
        if not go_sources_arg:
            assert self.go2obj_orig, \
                "go2obj MUST BE PRESENT IF go_sources IS NOT: {O}".format(O=self.go2obj_orig)
            self.go_sources = set(self.go2obj_orig)
            self.go2obj = self.go2obj_orig
            sys.stdout.write("**NOTE: {N:,} SOURCE GO IDS\n".format(N=len(self.go_sources)))
            return
        # GO sources provided
        go_sources = self._init_go_sources(go_sources_arg, self.go2obj_orig)
        # Create new go2obj_user subset matching GO sources
        # Fill with source and parent GO IDs and alternate GO IDs
        go2obj_user = {}
        objrel = CurNHigher(relationships_arg, self.go2obj_orig)
        objrel.get_id2obj_cur_n_high(go2obj_user, go_sources)
        # Add additional GOTerm information, if needed for user task
        kws_gos = {k:v for k, v in self.kws.items() if k in self.kws_aux_gos}
        if kws_gos:
            self._add_goterms_kws(go2obj_user, kws_gos)
        self.go_sources = go_sources
        self.go2obj = go2obj_user

    def _add_goterms_kws(self, go2obj_user, kws_gos):
        """Add more GOTerms to go2obj_user, if requested and relevant."""
        if 'go2color' in kws_gos:
            go2color = kws_gos['go2color']
            if go2color is not None:
                for goid in go2color.keys():
                    self._add_goterms(go2obj_user, goid)

    def _add_goterms(self, go2obj_user, goid):
        """Add alt GO IDs to go2obj subset, if requested and relevant."""
        goterm = self.go2obj_orig[goid]
        if goid != goterm.id and goterm.id in go2obj_user and goid not in go2obj_user:
            go2obj_user[goid] = goterm

    def _init_go_sources(self, go_srcs, go2obj_arg):
        """Return GO sources which are present in GODag."""
        gos_user = set(go_srcs) if not isinstance(go_srcs, str) else set([go_srcs,])
        if 'children' in self.kws and self.kws['children']:
            gos_user |= get_leaf_children(gos_user, go2obj_arg)
        gos_godag = set(go2obj_arg)
        gos_source = gos_user.intersection(gos_godag)
        gos_missing = gos_user.difference(gos_godag)
        if not gos_missing:
            return gos_source
        sys.stdout.write("{N} GO IDs NOT FOUND IN GO DAG: {GOs}\n".format(
            N=len(gos_missing), GOs=" ".join([str(e) for e in gos_missing])))
        return gos_source


class InitFields:
    """Initialize print attributes and namedtuple fields."""

    exp_keys = set(['rcntobj', 'tcntobj', 'go2nt', 'go2letter'])

    def __init__(self, ini_main, **kws):
        self.go2obj = ini_main.go2obj
        self.kws = get_kwargs(kws, self.exp_keys, None)
        if 'rcntobj' not in kws:
            self.kws['rcntobj'] = True
        self.kw_elems = self._init_kwelems()
        self.relationships = ini_main.relationships
        self.prt_flds = self._init_prt_flds()

    def get_rcntobj(self):
        """Return None or user-provided CountRelatives object."""
        # rcntobj value in kws can be: None, False, True, CountRelatives object
        if 'rcntobj' in self.kws:
            rcntobj = self.kws['rcntobj']
            if isinstance(rcntobj, CountRelatives):
                return rcntobj
            return CountRelatives(
                self.go2obj,  # Subset go2obj contains only items needed by go_sources
                self.relationships,
                dcnt='dcnt' in self.kw_elems,
                go2letter=self.kws.get('go2letter'))
        return None

    def get_go2nt_all(self, rcntobj):
        """For each GO id, put all printable fields in one namedtuple."""
        if 'go2nt' in self.kws:
            go2nt = self.kws['go2nt']
            return {go:go2nt[go] for go in self.go2obj}
        return self._get_go2nt_all(rcntobj)

    def _init_prt_flds(self):
        """Return the print fields in the go2nt namedtuple."""
        # Create namedtuple fields or copy namedtuple fields
        if 'go2nt' not in self.kws:
            return self.__init_prt_flds()
        return next(iter(self.kws['go2nt'].values()))._asdict()

    def __init_prt_flds(self):
        """Return the print fields in the go2nt namedtuple."""
        prt_flds = ['NS', 'level', 'depth']
        if self.relationships:
            prt_flds.append('reldepth')
        prt_flds.extend(['GO', 'alt', 'GO_name'])
        if 'dcnt' in self.kw_elems:
            prt_flds.append('dcnt')
        if 'D1' in self.kw_elems:
            prt_flds.append('D1')
        if 'tcnt' in self.kw_elems:
            prt_flds.append('tcnt')
            prt_flds.append('tfreq')
            prt_flds.append('tinfo')
        if self.relationships:
            prt_flds.append('childcnt')
            prt_flds.append('REL')
            prt_flds.append('REL_short')
            prt_flds.append('rel')
        prt_flds.append('id')
        return prt_flds

    def get_prt_fmt(self, alt=False):
        """Return the format for printing GO named tuples and their related information."""
        # prt_fmt = [ #                                                        rcnt
        #     '{GO} # {NS}  L{level:02} D{depth:02} {GO_name}',
        #     '{GO} # {NS} {dcnt:6,} L{level:02} D{depth:02} {D1:5} {GO_name}']
        prt_fmt = []
        if alt:
            prt_fmt.append('{GO}{alt:1}')
        else:
            prt_fmt.append('{GO}')
        prt_fmt.append('# {NS}')
        if 'dcnt' in self.prt_flds:
            prt_fmt.append('{dcnt:5}')
        if 'childcnt' in self.prt_flds:
            prt_fmt.append('{childcnt:3}')
        if 'tcnt' in self.prt_flds:
            prt_fmt.append("{tcnt:7,}")
        if 'tfreq' in self.prt_flds:
            prt_fmt.append("{tfreq:8.6f}")
        if 'tinfo' in self.prt_flds:
            prt_fmt.append("{tinfo:5.2f}")
        prt_fmt.append('L{level:02} D{depth:02}')
        if self.relationships:
            prt_fmt.append('R{reldepth:02}')
        if 'D1' in self.prt_flds:
            prt_fmt.append('{D1:5}')
        if 'REL' in self.prt_flds:
            prt_fmt.append('{REL}')
            prt_fmt.append('{rel}')
        prt_fmt.append('{GO_name}')
        return " ".join(prt_fmt)

    def _get_go2nt_all(self, rcntobj):
        """For each GO id, put all printable fields in one namedtuple."""
        ### tic = timeit.default_timer()
        go2nt = {}
        ntobj = cx.namedtuple("NtGo", " ".join(self.prt_flds))
        ### tic = _rpt_hms(tic, "GoSubDag: _Init::get_go2nt")
        tcntobj = self.kws['tcntobj'] if 'tcntobj' in self.kws else None
        b_tcnt = tcntobj is not None
        # b_rcnt = rcntobj is not None and rcntobj
        objrelstr = RelationshipStr(self.relationships)
        for goid, goobj in self.go2obj.items():
            ns_go = NAMESPACE2NS[goobj.namespace]
            fld2vals = {
                'NS' : ns_go,
                'level' : goobj.level,
                'depth' : goobj.depth,
                'GO' : goid,
                'alt' : '' if goid == goobj.id else 'a',
                'id' : goobj.id,
                'GO_name' : goobj.name}
            if 'dcnt' in self.kw_elems:
                fld2vals['dcnt'] = rcntobj.go2dcnt[goid]
            if 'D1' in self.kw_elems:
                fld2vals['D1'] = rcntobj.get_d1str(goobj)
            if b_tcnt:
                tcnt = tcntobj.gocnts.get(goid, 0)
                num_ns = float(tcntobj.aspect_counts[goobj.namespace])
                tfreq = float(tcnt)/num_ns if num_ns != 0 else 0
                fld2vals['tcnt'] = tcnt
                fld2vals['tfreq'] = tfreq
                fld2vals['tinfo'] = 0.0 - math.log(tfreq) if tfreq else 0
            if self.relationships:
                fld2vals['childcnt'] = len(goobj.children)
                fld2vals['reldepth'] = goobj.reldepth
                fld2vals['REL'] = objrelstr.str_relationships(goobj)
                fld2vals['REL_short'] = objrelstr.str_rel_short(goobj)
                fld2vals['rel'] = objrelstr.str_relationships_rev(goobj)
            go2nt[goid] = ntobj(**fld2vals)
        ### tic = _rpt_hms(tic, "GoSubDag: _Init::get_go2nt")
        return go2nt

    def _init_kwelems(self):
        """Init set elements."""
        ret = set()
        if 'rcntobj' in self.kws:
            ret.add('dcnt')
            ret.add('D1')
        if 'tcntobj' in self.kws:
            ret.add('tcnt')
            ret.add('tfreq')
            ret.add('tinfo')
        return ret


# Copyright (C) 2016-2020, DV Klopfenstein, H Tang, All rights reserved.
