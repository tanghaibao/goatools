"""Test the loading of the optional GO term fields."""

__copyright__ = "Copyright (C) 2010-2018, DV Klopfenstein, H Tang, All rights reserved."
__author__ = "DV Klopfenstein"


import os
import sys
import re
import collections as cx
import timeit
import datetime
from goatools.obo_parser import GODag
from goatools.obo_parser import OBOReader
from goatools.base import download_go_basic_obo


class OptionalAttrs(object):
    """Holds data for GO relationship test."""

    repo = os.path.join(os.path.dirname(os.path.abspath(__file__)), "../..")
    cmpfld = re.compile(r'^(\S+): (\S.*\S)\s*$')  # Field line pattern

    def __init__(self, fin_obo, opt_field=None):
        self.opt = opt_field  # None causes all fields to read to exp dict
        self.obo = os.path.join(self.repo, fin_obo)
        self._init_dnld_dag()
        self.go2obj = {o.id:o for o in self.load_dag(opt_field).values()}
        self.dcts = self._init_go2dct()  # dagdct go2dct typdefdct flds
        self.go2dct = {go:d for go, d in self.dcts['go2dct'].items() if go in self.go2obj}
        self.num_tot = len(self.go2obj)

    def prt_summary(self, prt=sys.stdout):
        """Print percentage of GO IDs that have a specific relationship."""
        prt.write("\nThese fields appear at least once on a GO term:\n")
        # Ex: 28,951 of 44,948 (64%) GO IDs has field(synonym)
        for relname, cnt in self._get_cnts_gte1().most_common():
            self._prt_perc(cnt, "relationship:{R}".format(R=relname), prt)
        prt.write("\nMaximum number of fields on a GO term:\n")
        for fld, maxqty in sorted(self._get_cnts_max().items(), key=lambda t: t[1]):
            prt.write("    {MAX:3} {MRK} {FLD}\n".format(
                MAX=maxqty, MRK=self._get_fldmrk(fld), FLD=fld))

    @staticmethod
    def _get_fldmrk(fld):
        """Get a mark for each field indicating if it is required or optional"""
        if fld in OBOReader.attrs_req or fld == 'def':
            return 'REQ'
        if fld in OBOReader.attrs_scalar:
            return 'str'
        if fld in OBOReader.attrs_set:
            return 'set'
        if fld == 'relationship':
            return 'rel'
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
        for goid, dct in self.go2dct.items():
            act_set = getattr(self.go2obj[goid], opt, None)
            if opt in dct:
                exp_num = len(dct[opt])
                act_num = len(act_set)
                assert exp_num == act_num, "SET EXP({}) ACT({}) {}:\nEXP({})\nACT({})".format(
                    exp_num, act_num, goid, dct[opt], act_set)
            else:
                assert act_set is None

    def chk_cnt_relationship(self, opt='relationship'):
        """Check that all GO IDs that should have relationships do have relationships."""
        for goid, dct in self.go2dct.items():
            act_rel2gos = getattr(self.go2obj[goid], opt, None)
            if opt in dct:
                exp_gos = set(s.split()[1] for s in dct[opt])
                act_gos = set(go for gos in act_rel2gos.values() for go in gos)
                assert exp_gos == act_gos, "EXP({}) ACT({}) {}:\nEXP({})\nACT({})".format(
                    len(exp_gos), len(act_gos), goid, exp_gos, act_gos)
            else:
                assert act_rel2gos is None

    def chk_cnt_go(self, opt):
        """Check that all GO IDs that should have relationships do have relationships."""
        if opt == 'def':
            opt = 'defn'
        num_exp = sum(1 for d in self.go2dct.values() if opt in d)
        go_act = set(o.id for o in self.go2obj.values() if hasattr(o, opt))
        assert num_exp == len(go_act), "GO CNTS({}): EXP({}) ACT({})".format(
            opt, num_exp, len(go_act))

    def _init_go2dct(self):
        """Get a list of GO IDs that have relationships."""
        dagdct = {}
        go2dct = {}
        typedefdct = {}
        flds = set()
        with open(self.obo) as ifstrm:
            rec = {}
            for line in ifstrm:
                line = line.rstrip()
                # Beginning of GO record
                if line[:7] == "id: GO:":
                    assert not rec, "NOW({}) WAS({})".format(line[4:], rec)
                    rec = {'GO':line[4:]}
                # End of GO record
                elif not line:
                    if rec:  # and option is None or option in rec:
                        # 'Definition' is specified in obo as 'def' and in Python by 'defn'
                        if 'def' in rec:
                            rec['defn'] = rec['def']
                        go2dct[rec['GO']] = rec
                    rec = {}
                # Middle of GO record
                elif rec:
                    # r'^(\S+): (\S.*\S)\s*$')  # Field line pattern
                    mtch = self.cmpfld.match(line)
                    if mtch:
                        fld = mtch.group(1)
                        flds.add(fld)
                        if fld not in rec:
                            rec[fld] = set()
                        rec[fld].add(mtch.group(2))
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
