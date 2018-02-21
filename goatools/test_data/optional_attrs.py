"""Test the loading of the optional GO term fields."""
# https://owlcollab.github.io/oboformat/doc/GO.format.obo-1_4.html

__copyright__ = "Copyright (C) 2010-2018, DV Klopfenstein, H Tang, All rights reserved."
__author__ = "DV Klopfenstein"


import os
import sys
import re
import collections as cx
import timeit
import datetime
from goatools.obo_parser import GODag
from goatools.base import download_go_basic_obo


class OptionalAttrs(object):
    """Holds data for GO relationship test."""

    repo = os.path.join(os.path.dirname(os.path.abspath(__file__)), "../..")
    cmpfld = re.compile(r'^(\S+)\s*:\s*(\S.*\S)\s*$')  # Field line pattern
    exp_scopes = set(['EXACT', 'BROAD', 'NARROW', 'RELATED'])
    exp_xrefpat = re.compile(r'^\S+:\S+$')
    # Required attributes are always loaded
    exp_req = set(['name', 'id', 'is_obsolete', 'namespace', 'alt_ids', 'is_a'])
    # Generated attributes
    exp_gen = set(['level', 'depth', 'parents', 'children', '_parents'])

    attrs_req = set(['id', 'namespace', 'name', 'is_a', 'alt_ids'])
    attrs_scalar = set(['id', 'namespace', 'name', 'def', 'comment'])
    attrs_set = set(['xref', 'subset', 'alt_id'])

    def __init__(self, fin_obo, opt_field=None):
        self.opt = opt_field  # None causes all fields to read to exp dict
        self.obo = os.path.join(self.repo, fin_obo)
        self._init_dnld_dag()
        self.go2obj = {o.id:o for o in self.load_dag(opt_field).values()}
        self.dcts = self._init_go2dct()  # dagdct go2dct typdefdct flds
        self.go2dct = {go:d for go, d in self.dcts['go2dct'].items() if go in self.go2obj}
        self.num_tot = len(self.go2obj)
        self._chk_required()
        self._chk_parents()
        self._set_exp_children()
        self._chk_children()

    def prt_summary(self, prt=sys.stdout):
        """Print percentage of GO IDs that have a specific relationship."""
        prt.write("\nThese fields appear at least once on a GO term:\n")
        # Ex: 28,951 of 44,948 (64%) GO IDs has field(synonym)
        for relname, cnt in self._get_cnts_gte1().most_common():
            self._prt_perc(cnt, relname, prt)
        prt.write("\nMaximum number of fields on a GO term:\n")
        for fld, maxqty in sorted(self._get_cnts_max().items(), key=lambda t: t[1]):
            prt.write("    {MAX:3} {MRK} {FLD}\n".format(
                MAX=maxqty, MRK=self._get_fldmrk(fld), FLD=fld))

    def _chk_parents(self):
        """Check that parent relationships."""
        for goobj in self.go2obj.values():
            exp_dct = self.go2dct[goobj.id]
            if 'is_a' in exp_dct:
                # pylint: disable=protected-access
                exp_parents = exp_dct['is_a']
                act_parents = goobj._parents
                assert exp_parents == act_parents
            else:
                assert not goobj.parents

    def _chk_children(self):
        """Check that parent relationships."""
        for goobj in self.go2obj.values():
            exp_dct = self.go2dct[goobj.id]
            if '_children' in exp_dct:
                exp_children = exp_dct['_children']
                act_children = set(o.id for o in goobj.children)
                assert exp_children == act_children
            else:
                assert not goobj.children

    def _set_exp_children(self):
        """Fill expected child GO IDs."""
        # Initialize empty sets for child GO IDs
        for exp_dct in self.go2dct.values():
            exp_dct['_children'] = set()
        # Loop thru all GO IDs
        for goid_child, exp_dct in self.go2dct.items():
            if 'is_a' in exp_dct:
                # Add current GO ID to all of it's parents' set of children
                for goid_parent in exp_dct['is_a']:
                    self.go2dct[goid_parent]['_children'].add(goid_child)

    def _chk_required(self):
        """Check the required attributes."""
        for goid, goobj in self.go2obj.items():
            godct = self.go2dct[goid]
            assert goobj.id == godct['GO']
            assert goobj.namespace == next(iter(godct['namespace'])), godct
            assert goobj.name == next(iter(godct['name']))
            self._chk_is_obsolete(goobj, godct)
            self._chk_alt_ids(goobj, godct)

    @staticmethod
    def _chk_alt_ids(goobj, godct):
        """Check 'alt_ids' required attribute."""
        if 'alt_id' in godct:
            assert godct['alt_id'] == goobj.alt_ids
        else:
            assert not goobj.alt_ids

    @staticmethod
    def _chk_is_obsolete(goobj, godct):
        """Check 'is_obsolete' required attribute."""
        act_obso = getattr(goobj, 'is_obsolete', None)
        if act_obso:
            assert 'is_obsolete' in godct, "EXP({})\nACT({})".format(
                godct, getattr(goobj, 'is_obsolete', None))
        else:
            assert 'is_obsolete' not in godct, "EXP({})\nACT({})".format(
                godct, getattr(goobj, 'is_obsolete', None))

    def chk_no_optattrs(self):
        """Check that only the optional attributes requested are the attributes implemented."""
        # name is_obsolete namespace id alt_ids
        # level namespace depth parents children _parents
        exp_flds = self.exp_req.union(self.exp_gen)
        for goobj in self.go2obj.values():
            assert not set(vars(goobj).keys()).difference(exp_flds)
            # print(vars(goobj).keys())
            # print(" ".join(vars(goobj).keys()))

    def chk_xref(self, prt=None):
        """Check synonyms."""
        # Get GO IDs which are expected to have synonyms
        goids = set(go for go, d in self.go2dct.items() if 'xref' in d)
        for goid in goids:
            goobj = self.go2obj[goid]
            xrefs = getattr(goobj, 'xref', None)
            assert xrefs is not None, "{GO} MISSING XREF".format(GO=goid)
            # Iterate through list of synonym data stored in named tuples
            for dbxref in xrefs:
                if prt is not None:
                    prt.write("{GO} {DBXREF}\n".format(GO=goid, DBXREF=dbxref))
                assert self.exp_xrefpat.match(dbxref), "INVALID XREF FORMAT"

    def chk_synonyms(self, prt=None):
        """Check synonyms."""
        # Get GO IDs which are expected to have synonyms
        goids = set(go for go, d in self.go2dct.items() if 'synonym' in d)
        for goid in goids:
            goobj = self.go2obj[goid]
            ntsyns = getattr(goobj, 'synonym', None)
            assert ntsyns is not None, "{GO} MISSING SYNONYM".format(GO=goid)
            # Iterate through list of synonym data stored in named tuples
            for ntsyn in ntsyns:
                if prt is not None:
                    prt.write("{GO} {NT}\n".format(GO=goid, NT=ntsyn))
                assert ntsyn.text, "SYNONYM CANNOT BE EMPTY"
                assert ntsyn.scope in self.exp_scopes, "INVALID SYNONYM SCOPE"
                for dbxref in ntsyn.dbxrefs:
                    assert self.exp_xrefpat.match(dbxref), "INVALID SYNONYM DBXREF"

    def _get_fldmrk(self, fld):
        """Get a mark for each field indicating if it is required or optional"""
        #pylint: disable=too-many-return-statements
        if fld in self.attrs_req:
            return 'REQ'
        if fld == 'def':
            return 'str'
        if fld in self.attrs_scalar:
            return 'str'
        if fld in self.attrs_set:
            return 'set'
        if fld == 'relationship':
            return 'rel'
        if fld == 'synonym':
            return 'syn'
        if fld == 'xref':
            return 'xrf'
        raise RuntimeError("UNEXPECTED FIELD({})".format(fld))

    def _prt_perc(self, num_rel, name, prt=sys.stdout):
        """Print percentage of GO IDs that have a specific relationship."""
        prt.write("    {N:6,} of {M:,} ({P:3.0f}%) GO IDs has field({A})\n".format(
            N=num_rel, M=self.num_tot, P=float(num_rel)/self.num_tot*100, A=name))

    def _get_cnts_max(self):
        """Get the maximum count of times a specific relationship was seen on a GO."""
        fld2qtys = cx.defaultdict(set)
        flds = self.dcts['flds']
        for recdct in self.go2dct.values():
            for opt in flds:
                if opt in recdct:
                    fld2qtys[opt].add(len(recdct[opt]))
        return {f:max(qtys) for f, qtys in fld2qtys.items()}

    def _get_cnts_gte1(self):
        """Get counts of if a specific relationship was seen on a GO."""
        ctr = cx.Counter()
        flds = self.dcts['flds']
        for recdct in self.go2dct.values():
            for opt in flds:
                if opt in recdct:
                    ctr[opt] += 1
        return ctr

    def chk_cnt_set(self, opt):
        """For each GO ID, check that actual count of a set attr equals expected count."""
        errpat = "SET EXP({EXP}) ACT({ACT}) {GO}\n{DESC}\n:\nEXP:\n{Es}\n\nACT:\n{As}"
        for goid, dct in self.go2dct.items():
            act_set = getattr(self.go2obj[goid], opt, None)
            if opt in dct:
                exp_num = len(dct[opt])
                act_num = len(act_set)
                assert exp_num == act_num, errpat.format(
                    EXP=exp_num, ACT=act_num, GO=goid,
                    DESC=str(self.go2obj[goid]),
                    Es="\n".join(sorted(dct[opt])),
                    As="\n".join(sorted(act_set)))
            else:
                assert act_set is None

    def chk_relationships(self):
        """Check that all GO IDs that should have relationships do have relationships."""
        opt = 'relationship'
        for goid, dct in self.go2dct.items():
            act_rel2goobjs = getattr(self.go2obj[goid], opt, None)
            if opt in dct:
                exp_goids = set(s.split()[1] for s in dct[opt])
                act_goids = set(goobj.id for goobjs in act_rel2goobjs.values() for goobj in goobjs)
                assert exp_goids == act_goids, "EXP({}) ACT({}) {}:\nEXP({})\nACT({})".format(
                    len(exp_goids), len(act_goids), goid, exp_goids, act_goids)
            else:
                assert act_rel2goobjs is None

    def chk_cnt_go(self, opt):
        """Check that all GO IDs that should have relationships do have relationships."""
        if opt == 'def':
            opt = 'defn'
        num_exp = sum(1 for d in self.go2dct.values() if opt in d)
        go_act = set(o.id for o in self.go2obj.values() if hasattr(o, opt))
        assert num_exp == len(go_act), "GO CNTS({}): EXP({}) ACT({})".format(
            opt, num_exp, len(go_act))

    def _init_go2dct(self):
        """Create a dict of GO fields for use as expected results during test."""
        dagdct = {}
        go2dct = {}
        typedefdct = {}
        flds = set()
        with open(self.obo) as ifstrm:
            rec = {}
            for line in ifstrm:
                line = line.rstrip()
                # End of GO record
                if not line:
                    if rec:  # and option is None or option in rec:
                        # 'Definition' is specified in obo as 'def' and in Python by 'defn'
                        if 'def' in rec:
                            rec['defn'] = rec['def']
                        go2dct[rec['GO']] = rec
                    rec = {}
                else:
                    mtch = self.cmpfld.match(line)
                    if mtch:
                        fld = mtch.group(1)
                        val = mtch.group(2)

                        # Beginning of GO record
                        if fld == "id":
                            assert not rec, "NOW({}) WAS({})".format(line, rec)
                            rec = {'GO':val}
                        # Middle of GO record
                        elif rec:
                            flds.add(fld)
                            if fld not in rec:
                                rec[fld] = set()
                            # Strip comment if it exists
                            loc = val.find(' ! ')
                            if loc != -1:
                                val = val[:loc]
                            # Add value
                            rec[fld].add(val)
        return {'dagdct':dagdct, 'go2dct':go2dct, 'typedefdct':typedefdct, 'flds':flds}

    def _init_dnld_dag(self):
        """If dag does not exist, download it."""
        if not os.path.exists(self.obo):
            download_go_basic_obo(self.obo, loading_bar=None)

    def load_dag(self, opt_fields=None):
        """Run numerous tests for various self.reports."""
        tic = timeit.default_timer()
        dag = GODag(self.obo, opt_fields)
        toc = timeit.default_timer()
        msg = "Elapsed HMS for OBO DAG load: {HMS} OPTIONAL_ATTR({O})\n".format(
            HMS=str(datetime.timedelta(seconds=(toc-tic))), O=opt_fields)
        sys.stdout.write(msg)
        return dag


# Copyright (C) 2010-2018, DV Klopfenstein, H Tang, All rights reserved.
