#!/usr/bin/env python
"""Test the loading of the optional GO term field, 'relationship'.

In go-basic.obo fmt(1.2) rel(2018-02-11) 15,329 out of 47,120 GO Terms have relationships:
   6904 part_of
   3226 regulates
   2802 negatively_regulates
   2784 positively_regulates


"""

import os
import sys
import collections as cx
import timeit
import datetime
from goatools.obo_parser import GODag
from goatools.base import download_go_basic_obo


# GO:0007608 L06 D06 sensory perception of smell
# GO:0050911 L05 D05 detection of chemical stimulus involved in sensory perception of smell
# GO:0007260

def test_defn():
    """Test loading optional GO term field, 'relationship'."""
    obj = _Run("go-basic.obo", "relationship")
    # obj.load_dag(obo)  # 0:00:03.017731 (HMS for GO-DAG load w/no optional attrs)
    # Check that all GO IDs that should have relationships do have relationships.
    obj.prt_summary()
    obj.chk_cnt_go()
    obj.chk_cnt_rel()


class _Run(object):
    """Holds data for GO relationship test."""

    repo = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..")

    def __init__(self, fin_obo, opt_field):
        self.opt = opt_field
        self.obo = os.path.join(self.repo, fin_obo)
        self.go2obj = self._init_go2obj(opt_field)
        self.num_tot = len(self.go2obj)
        self.go2relset_exp = self._init_go2relset_exp()  # For EXPECTED results

    def prt_summary(self, prt=sys.stdout):
        """Print percentage of GO IDs that have a specific relationship."""
        prt.write("\nSummary of frequency of relationship on GO terms:\n")
        num_rel = len(self.go2relset_exp)
        self._prt_perc(num_rel, self.opt, prt)
        for relname, cnt in self._get_cnts().most_common():
            self._prt_perc(cnt, "relationship:{R}".format(R=relname), prt)

    def _prt_perc(self, num_rel, name, prt=sys.stdout):
        """Print percentage of GO IDs that have a specific relationship."""
        prt.write("    {N:6,} of {M:,} ({P:2.0f}%) GO IDs has field({A})\n".format(
            N=num_rel, M=self.num_tot, P=float(num_rel)/self.num_tot*100, A=name))

    def _get_cnts(self):
        """Get counts of if a specific relationship was seen on a GO."""
        ctr = cx.Counter()
        for goid in self.go2relset_exp:
            rel2gos = getattr(self.go2obj[goid], self.opt)
            # num_gos = len(set(go for gos in rel2gos.values() for go in gos))
            # if num_gos > 1:
            #     print("RRRRRRRRRRRRRRR", goid, rel2gos)
            for rel in rel2gos:
                ctr[rel] += 1 
        return ctr

    def chk_cnt_rel(self):
        """Check that all GO IDs that should have relationships do have relationships."""
        for goid, relstrs in self.go2relset_exp.items():
            exp_num = len(relstrs)
            act_rel2gos = self.go2obj[goid].relationship
            act_num = sum(len(gos) for gos in act_rel2gos.values())
            assert exp_num == act_num, "EXP({}) ACT({}) {}:\nEXP({})\nACT({})".format(
                exp_num, act_num, goid, relstrs, act_rel2gos)

    def chk_cnt_go(self):
        """Check that all GO IDs that should have relationships do have relationships."""
        num_exp = len(self.go2relset_exp)
        go_act = set(o.id for o in self.go2obj.values() if hasattr(o, 'relationship'))
        assert num_exp == len(go_act), "GO CNTS: EXP({}) ACT({})".format(
            num_exp, len(go_act))

    def _init_go2relset_exp(self):
        """Get a list of GO IDs that have relationships."""
        go2relset = {}
        with open(self.obo) as ifstrm:
            rec = {}
            for line in ifstrm:
                line = line.rstrip()
                if line[:7] == "id: GO:":
                    assert not rec, "NOW({}) WAS({})".format(line[4:], rec)
                    rec = {'GO':line[4:]}
                elif line[:14] == "relationship: ":
                    if 'rel' not in rec:
                        rec['rel'] = set()
                    rec['rel'].add(line[14:])
                elif not line:
                    if rec and 'rel' in rec:
                        go2relset[rec['GO']] = rec['rel']
                    rec = {}
        return go2relset

    def _init_go2obj(self, opt_field):
        self._init_dnld_dag()
        # _load_dag(obo)  # 0:00:03.017731
        return {o.id:o for o in self.load_dag(opt_field).values()}

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


if __name__ == '__main__':
    test_defn()
