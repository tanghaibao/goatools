"""Test the loading of the optional GO term fields."""
# https://owlcollab.github.io/oboformat/doc/GO.format.obo-1_4.html

__copyright__ = "Copyright (C) 2010-2018, DV Klopfenstein, H Tang, All rights reserved."
__author__ = "DV Klopfenstein"


import os
import sys
import re
import timeit
import collections as cx
from goatools.godag.prttime import GoDagTimed
from goatools.godag.prttime import prt_hms


class OptionalAttrs(object):
    """Holds data for GO relationship test."""

    repo = os.path.join(os.path.dirname(os.path.abspath(__file__)), "../..")
    cmpfld = re.compile(r'^(\S+)\s*:\s*(\S.*\S)\s*$')  # Field line pattern
    exp_scopes = set(['EXACT', 'BROAD', 'NARROW', 'RELATED'])
    exp_xrefpat = re.compile(r'^\S+:\S+$')
    # Required attributes are always loaded
    exp_req = set(['name', 'item_id', 'is_obsolete', 'namespace', 'alt_id', 'is_a', 'is_obsolete'])
    # Generated attributes
    exp_gen = set(['level', 'depth', 'parents', 'children', '_parents'])
    exp_relationships = set(['part_of',
                             'regulates', 'negatively_regulates', 'positively_regulates'])

    attrs_scalar = set(['item_id', 'namespace', 'name', 'def', 'comment'])
    attrs_set = set(['xref', 'subset', 'alt_id'])

    def __init__(self, fin_obo, opt_field=None, keep_alt_ids=False):
        self.opt = opt_field  # None causes all fields to read to exp dict
        self.obo = os.path.join(self.repo, fin_obo)
        self.go2obj = GoDagTimed(self.obo, opt_field, keep_alt_ids).go2obj
        self.dcts = self._init_go2dct()  # go2dct typdefdct flds
        self.go2dct = {go:d for go, d in self.dcts['go2dct'].items() if go in self.go2obj}
        self.num_tot = len(self.go2obj)
        self._chk_required()
        self._chk_parents()
        self._set_exp_children()
        self._chk_children()

    def chk_get_goterms_upper(self):
        """Check that GOTerm's 'get_upper' returns parents and relationships."""
        tic = timeit.default_timer()
        for goterm in self.go2obj.values():
            goids_act = set(o.item_id for o in goterm.get_goterms_upper())
            goids_exp = self._get_goterms_upper(goterm.item_id)
            assert goids_act == goids_exp
        prt_hms(tic, "get_goterms_upper")

    def chk_get_goterms_lower(self):
        """Check that GOTerm's 'get_lower' returns parents and relationships."""
        tic = timeit.default_timer()
        for goterm in self.go2obj.values():
            goids_act = set(o.item_id for o in goterm.get_goterms_lower())
            goids_exp = self._get_goterms_lower(goterm.item_id)
            assert goids_act == goids_exp, "{GO} EXP({E}) ACT({A})".format(
                GO=goterm.item_id, E=goids_exp, A=goids_act)
        prt_hms(tic, "get_goterms_lower")

    def _get_goterms_upper(self, goid):
        """Get expected GO IDs returned by GOTerm's 'get_goterms_upper'."""
        goids_exp = set()
        dct = self.go2dct[goid]
        if 'is_a' in dct:
            goids_exp.update(dct['is_a'])
        if 'relationship' in dct:
            for rel_go in dct['relationship']:
                goids_exp.add(rel_go.split()[1])
        return goids_exp

    def _get_goterms_lower(self, goid):
        """Get expected GO IDs returned by GOTerm's 'get_goterms_lower'."""
        goids_exp = set()
        dct = self.go2dct[goid]
        if 'is_a_rev' in dct:
            goids_exp.update(dct['is_a_rev'])
        if 'relationship_rev' in dct:
            for rel_gos in dct['relationship_rev'].values():
                goids_exp.update(rel_gos)
        return goids_exp

    def chk_relationships_rev(self, reltype='part_of', prt=None):
        """Check reciprocal relationships. Print all GO pairs in one type of relationship."""
        spc = " "*len(reltype)
        rec2revs = cx.defaultdict(set)
        for rec in sorted(self.go2obj.values(), key=lambda o: o.namespace):
            reldct = rec.relationship
            if reltype in reldct:
                if prt is not None:

                    prt.write("{SPC} {GO}\n".format(SPC=spc, GO=str(rec)))
                for related_to in reldct[reltype]:
                    rec2revs[related_to].add(rec)
                    if prt is not None:
                        prt.write("{RELTYPE} {GO}\n".format(RELTYPE=reltype, GO=str(related_to)))
                if prt is not None:
                    prt.write("\n")
        for rec, exp_revs in sorted(rec2revs.items(), key=lambda t: t[0].namespace):
            if prt is not None:
                prt.write("    {SPC} {GO}\n".format(SPC=spc, GO=str(rec)))
            assert rec.relationship_rev[reltype] == exp_revs
            for related_from in rec.relationship_rev[reltype]:
                if prt is not None:
                    prt.write("rev {RELTYPE} {GO}\n".format(RELTYPE=reltype, GO=str(related_from)))
            if prt is not None:
                prt.write("\n")

    def chk_str(self, attr):
        """Check that expected scalar value matches actual string value."""
        for goid, rec in self.go2obj.items():
            # A string data member must always be present, even if the value is ""
            act_str = getattr(rec, attr)
            exp_dct = self.go2dct[goid]
            # Expected string equals actual string?
            if attr in exp_dct:
                exp_str = next(iter(exp_dct[attr]))
                assert exp_str == act_str, "{} EXP({}) ACT({})".format(
                    goid, exp_str, act_str)
            # If there is no expected string, is actual string ""?
            else:
                assert act_str == ""

    def prt_summary(self, prt=sys.stdout):
        """Print percentage of GO IDs that have a specific relationship."""
        sep = "\n-----------------------------------------------------------\n"
        flds_seen = self.dcts['flds']
        fld_cnts_go = self._get_cnts_gte1(self.go2dct.values())
        prt.write("{SEP}GO TERM REQUIRED FIELDS:\n".format(SEP=sep))
        self._prt_summary(prt, fld_cnts_go, self.exp_req, self.go2dct.values())
        flds_seen = flds_seen.difference(self.exp_req)
        prt.write("{SEP}GO TERM OPTIONAL FIELDS:\n".format(SEP=sep))
        self._prt_summary(prt, fld_cnts_go, flds_seen, self.go2dct.values())
        flds_seen = flds_seen.difference(fld_cnts_go.keys())
        prt.write("{SEP}Typedef FIELDS:\n".format(SEP=sep))
        fld_cnts_typedef = self._get_cnts_gte1(self.dcts['typedefdct'].values())
        self._prt_summary(prt, fld_cnts_typedef, flds_seen, self.dcts['typedefdct'])
        flds_seen = flds_seen.difference(fld_cnts_typedef.keys())
        assert flds_seen == set(['consider', 'replaced_by']), "UNEXPECTED FIELDS({})".format(
            flds_seen)

    def _prt_summary(self, prt, fld_cnts, prt_flds, dcts):
        prt.write("\n    These fields appear at least once\n")
        # Ex: 28,951 of 44,948 (64%) GO IDs has field(synonym)
        for relname, cnt in fld_cnts.most_common():
            if prt_flds is None or relname in prt_flds:
                self._prt_perc(cnt, relname, len(dcts), prt)
        prt.write("\n    Maximum number of fields:\n")
        for fld, maxqty in sorted(self._get_cnts_max(dcts).items(), key=lambda t: t[1]):
            if prt_flds is None or fld in prt_flds:
                prt.write("        {MAX:3} {MRK} {FLD}\n".format(
                    MAX=maxqty, MRK=self._get_fldmrk(fld), FLD=fld))

    def _chk_parents(self):
        """Check parents."""
        for goobj in self.go2obj.values():
            exp_dct = self.go2dct[goobj.item_id]
            if 'is_a' in exp_dct:
                # pylint: disable=protected-access
                exp_parents = exp_dct['is_a']
                act_parents = goobj._parents
                assert exp_parents == act_parents
            else:
                assert not goobj.parents

    def _chk_children(self):
        """Check children."""
        for goobj in self.go2obj.values():
            exp_dct = self.go2dct[goobj.item_id]
            if '_children' in exp_dct:
                exp_children = exp_dct['_children']
                act_children = set(o.item_id for o in goobj.children)
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
            assert goobj.item_id == godct['GO']
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
        # name is_obsolete namespace item_id alt_ids
        # level namespace depth parents children _parents
        exp_flds = self.exp_req.union(self.exp_gen)
        obj1_exp0 = set(['id', 'alt_ids'])
        for goobj in self.go2obj.values():
            attrs = set(vars(goobj).keys()).difference(exp_flds)
            assert attrs == obj1_exp0, attrs

    def chk_xref(self, prt=None):
        """Check xrefs."""
        # Get GO IDs which are expected to have xrefs
        goids = set(go for go, d in self.go2dct.items() if 'xref' in d)
        for goid in goids:
            goobj = self.go2obj[goid]
            xrefs = getattr(goobj, 'xref', None)
            assert xrefs is not None, "{GO} MISSING XREF".format(GO=goid)
            # Iterate through list of xref data stored in named tuples
            for dbxref in xrefs:
                if prt is not None:
                    prt.write("{GO} {DBXREF}\n".format(GO=goid, DBXREF=dbxref))
                assert self.exp_xrefpat.match(dbxref), "INVALID XREF FORMAT({X})".format(
                    X=dbxref)

    def chk_synonyms(self, prt=None):
        """Check synonyms

        Example synonym and its storage in a namedtuple:
        synonym: "The other white meat" EXACT MARKETING_SLOGAN [MEAT:00324, BACONBASE:03021]
          text:     "The other white meat"
          scope:    EXACT
          typename: MARKETING_SLOGAN
          dbxrefs:  set(["MEAT:00324", "BACONBASE:03021"])
        """
        # Get GO IDs which are expected to have synonyms
        badnts = []
        for goid, dct_exp in self.go2dct.items():
            goobj = self.go2obj[goid]
            if 'synonym' in dct_exp:
                ntsyns = getattr(goobj, 'synonym', None)
                assert ntsyns is not None, "{GO} MISSING SYNONYM".format(GO=goid)
                # Iterate through list of synonym data stored in named tuples
                for ntsyn in ntsyns:
                    if prt is not None:
                        prt.write("{GO} {NT}\n".format(GO=goid, NT=ntsyn))
                    # Example:
                    assert ntsyn.text, "SYNONYM CANNOT BE EMPTY"
                    assert ntsyn.scope in self.exp_scopes, "INVALID SYNONYM SCOPE"
                    for dbxref in ntsyn.dbxrefs:
                        if not self.exp_xrefpat.match(dbxref):
                            badnts.append((goid, ntsyn))
                            print("**WARNING: INVALID FORMAT: DBXREF({D}) ON {GO}".format(
                                D=dbxref, GO=goid))
            else:
                assert goobj.synonym == []
        return badnts

    def _get_fldmrk(self, fld):
        """Get a mark for each field indicating if it is required or optional"""
        #pylint: disable=too-many-return-statements
        if fld in self.exp_req:
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

    @staticmethod
    def _prt_perc(num_rel, name, num_tot, prt=sys.stdout):
        """Print percentage of GO IDs that have a specific relationship."""
        prt.write("        {N:6,} of {M:,} ({P:3.0f}%) GO IDs has field({A})\n".format(
            N=num_rel, M=num_tot, P=float(num_rel)/num_tot*100, A=name))

    def _get_cnts_max(self, dcts):
        """Get the maximum count of times a specific relationship was seen on a GO."""
        fld2qtys = cx.defaultdict(set)
        flds = self.dcts['flds']
        for recdct in dcts:
            for opt in flds:
                if opt in recdct:
                    fld2qtys[opt].add(len(recdct[opt]))
        return {f:max(qtys) for f, qtys in fld2qtys.items()}

    def _get_cnts_gte1(self, record_dicts):
        """Get counts of if a specific relationship was seen on a GO."""
        ctr = cx.Counter()
        flds = self.dcts['flds']
        for recdct in record_dicts:
            for opt in flds:
                if opt in recdct:
                    ctr[opt] += 1
        return ctr

    def chk_set(self, opt):
        """Check that actual set contents match expected set contents."""
        errpat = "SET EXP({EXP}) ACT({ACT}) {GO}\n{DESC}:\nEXP:\n{Es}\n\nACT:\n{As}"
        for goid, dct in self.go2dct.items():
            act_set = getattr(self.go2obj[goid], opt, None)
            if opt in dct:
                exp_set = dct[opt]
                assert exp_set == act_set, errpat.format(
                    EXP=len(exp_set), ACT=len(act_set), GO=goid,
                    DESC=str(self.go2obj[goid].name),
                    Es="\n".join(sorted(exp_set)),
                    As="\n".join(sorted(act_set)))
            else:
                assert act_set == set(), "EXPECTED EMPTY SET FOR {O}: ACT({A})\n".format(
                    O=opt, A=act_set)

    def chk_relationships(self):
        """Expected relationship GO IDs should match actual relationship GO IDs."""
        for goid, dct in self.go2dct.items():
            act_rel2recs = getattr(self.go2obj[goid], 'relationship', None)
            if 'relationship' in dct:
                rel2gos = self._mk_exp_relatinship_sets(dct['relationship'])
                # Check if expected relationships and actual relationships are the same
                assert set(act_rel2recs.keys()) == set(rel2gos.keys()), "EXP({}) != ACT({})".format(
                    set(act_rel2recs.keys()), set(rel2gos.keys()))
                for rel, exp_goids in rel2gos.items():
                    # Expected relationships store GO IDs.
                    # Actual relationships store GO Terms.
                    act_goids = set(o.item_id for o in act_rel2recs[rel])
                    assert exp_goids == act_goids, "EXP({}) ACT({}) {}:\nEXP({})\nACT({})".format(
                        len(exp_goids), len(act_goids), goid, exp_goids, act_goids)
            else:
                assert act_rel2recs == {}, act_rel2recs

    def _mk_exp_relatinship_sets(self, relationship_str_set):
        """Transform a set of relationship strings into a dict of sets containing GO IDs."""
        rel2gos = cx.defaultdict(set)
        for rel_str in relationship_str_set:
            rel, goid = rel_str.split()
            assert rel in self.exp_relationships
            assert goid[:3] == "GO:" and goid[3:].isdigit()
            rel2gos[rel].add(goid)
        return rel2gos

    @staticmethod
    def add_is_a_rev(go2dct):
        """If there 'is_a' exists, add 'is_a_rev'."""
        for go_src, dct in go2dct.items():
            if 'is_a' in dct:
                for go_parent in dct['is_a']:
                    if 'is_a_rev' not in go2dct[go_parent]:
                        go2dct[go_parent]['is_a_rev'] = set()
                    go2dct[go_parent]['is_a_rev'].add(go_src)

    @staticmethod
    def add_relationship_rev(go2dct):
        """If there is a relationship, add 'relationship_rev'."""
        for go_src, dct in go2dct.items():
            if 'relationship' in dct:
                for rel in dct['relationship']:
                    reltype, go_dst = rel.split()
                    # print("RRRRRRRRR", go_src, reltype, go_dst)
                    if 'relationship_rev' not in go2dct[go_dst]:
                        go2dct[go_dst]['relationship_rev'] = {}
                    if reltype not in go2dct[go_dst]['relationship_rev']:
                        go2dct[go_dst]['relationship_rev'][reltype] = set()
                    go2dct[go_dst]['relationship_rev'][reltype].add(go_src)

    # pylint: disable=too-many-branches
    def _init_go2dct(self):
        """Create EXPECTED RESULTS stored in a dict of GO fields."""
        go2dct = {}
        # pylint: disable=unsubscriptable-object
        typedefdct = {}
        flds = set()
        with open(self.obo) as ifstrm:
            rec = {}
            rec_typedef = None
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
                    if rec_typedef is not None:
                        # Example rec_typedef:
                        #     {'xref': 'RO:0002212',
                        #      'name': 'negatively regulates',
                        #      'namespace': 'external',
                        #      'transitive_over': 'part_of',
                        #      'is_a': 'regulates',
                        #      'id': 'negatively_regulates'}
                        typedefdct[rec_typedef['item_id']] = rec_typedef
                        rec_typedef = None
                elif line[:9] == "[Typedef]":
                    rec_typedef = {}
                else:
                    mtch = self.cmpfld.match(line)
                    if mtch:
                        fld = mtch.group(1)
                        val = mtch.group(2)

                        # Beginning of GO record
                        if fld == "id":
                            assert not rec, "NOW({}) WAS({})".format(line, rec)
                            rec = {'GO':val, 'item_id':val}
                            flds.add('item_id')
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

                        if rec_typedef is not None:
                            if fld == 'id':
                                fld = 'item_id'
                            rec_typedef[fld] = val

        for dct in go2dct.values():
            if 'def' in dct:
                dct['defn'] = dct['def']
        self.add_relationship_rev(go2dct)
        self.add_is_a_rev(go2dct)
        return {'go2dct':go2dct, 'typedefdct':typedefdct, 'flds':flds}


# Copyright (C) 2010-2018, DV Klopfenstein, H Tang, All rights reserved.
