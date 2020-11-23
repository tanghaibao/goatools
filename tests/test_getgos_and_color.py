#!/usr/bin/env python
"""Test relationships ane depth counted w/relationshipd (reldepths) and plotting"""

__copyright__ = "Copyright (C) 2016-present, DV Klopfenstein, H Tang, All rights reserved."

from os.path import join
from os.path import exists
from collections import defaultdict
from goatools.base import get_godag
from goatools.cli.gosubdag_plot import PlotCli
from goatools.cli.gos_get import get_go2color
from goatools.godag.consts import RELATIONSHIP_SET
from tests.utils import REPO

# pylint: disable=line-too-long
GODAG = get_godag(join(REPO, "tests/data/i126/viral_gene_silence.obo"), optional_attrs={'relationship'})

def test_getgo_color_for_plots():
    """Test relationships ane depth counted w/relationshipd (reldepths) and plotting"""
    oplt = PlotCli(gosubdag=None, use_doc=False)

    kws_base = {
        'GO': {'GO:0060150'},
        'obo': join(REPO, 'tests/data/i126/viral_gene_silence.obo'),
        'go_color_file': join(REPO, 'tests/data/i126/viral_gene_silence.txt'),
    }
    okws = _Run(kws_base)

    okws.chk_r0(oplt.plot(GODAG, okws.merge(outfile='viral_r0.png')))

    oplt.plot(GODAG, okws.merge(
        outfile='viral_r1.png',
        relationships=RELATIONSHIP_SET))

    # 24 orange GO ancestors up part_of relationship
    # 12 green  GO ancestors up is_a    relationship
    # GO:0008150   biological_process
    # -GO:0051704  multi-organism process
    # --GO:0051707 response to other organism
    # --GO:0044419 interspecies interaction between organisms
    okws.chk_partof(oplt.plot(GODAG, okws.merge(outfile='viral_r_partof.png', relationships='part_of')))

    oplt.get_gosubdagplot(GODAG, okws.merge(relationships={'negatively_regulates',}))
    okws.chk_regp(oplt.plot(GODAG, okws.merge(outfile='viral_r_regp.png', relationships='positively_regulates')))

    okws.chk_reg(oplt.plot(GODAG, okws.merge(
        outfile='viral_r_regrp.png',
        relationships='regulates,positively_regulates')))

    okws.chk_reg(oplt.plot(GODAG, okws.merge(
        outfile='viral_r_regrn.png',
        relationships='regulates,negatively_regulates')))

    _, oplt1 = oplt.plot(GODAG, okws.merge(
        outfile='viral_r_regpn.png',
        relationships='positively_regulates,negatively_regulates'))

    _, oplt2 = oplt.plot(GODAG, okws.merge(
        outfile='viral_r_regrpn.png',
        relationships='regulates,positively_regulates,negatively_regulates'))

    assert okws._get_reldepth(oplt1.gosubdag) == 10
    assert okws._get_reldepth(oplt2.gosubdag) == 11
    assert set(oplt1.gosubdag.go2nt.keys()) == set(oplt2.gosubdag.go2nt.keys())


class _Run:
    """Merge kw dicts"""

    hex2color = {
        '#e6fad2': 'green',
        '#ffe5b4': 'orange',
        '#d2d2fa': 'purple',
        '#d2fafa': 'blue',
        '#fad2fa': 'red',
    }

    def __init__(self, kws_base):
        self.kws_base = kws_base
        self.goid = next(iter(kws_base['GO']))
        self.go2color = get_go2color(kws_base['go_color_file'])
        self.color2goids = self._init_color2goids()

    def merge(self, **kws_usr):
        """Merge kws with base kws"""
        kws = dict(self.kws_base)
        for key, val in kws_usr.items():
            if key == 'outfile':
                val = join(REPO, val)
            kws[key] = val
        return kws

    def chk_reg(self, arg):
        """Check relationship, part_of"""
        png, gosubdagplot = arg
        assert exists(png)
        go2p = gosubdagplot.gosubdag.rcntobj.go2ancestors
        exp_all = set.union(*[self.color2goids[c] for c in ['green', 'purple', 'red', 'blue']])
        self._cmp(go2p[self.goid], exp_all)
        assert self._get_reldepth(gosubdagplot.gosubdag) == 11

    def chk_regp(self, arg):
        """Check relationship, part_of"""
        png, gosubdagplot = arg
        assert exists(png)
        go2p = gosubdagplot.gosubdag.rcntobj.go2ancestors
        exp_all = set.union(*[self.color2goids[c] for c in ['green', 'purple', 'red']]).difference({'GO:0071704'})
        self._cmp(go2p[self.goid], exp_all)
        assert self._get_reldepth(gosubdagplot.gosubdag) == 10

    def chk_partof(self, arg):
        """Check relationship, part_of"""
        png, gosubdagplot = arg
        assert exists(png)
        go2p = gosubdagplot.gosubdag.rcntobj.go2ancestors
        exp_all = self.color2goids['green'].union(self.color2goids['orange'])
        self._cmp(go2p[self.goid], exp_all)
        assert self._get_reldepth(gosubdagplot.gosubdag) == 9
        
    @staticmethod
    def _get_source_term(gosubdag):
        """Get the source term in the GoSubDag"""
        goid = next(iter(gosubdag.go_sources))
        return gosubdag.go2obj[goid]

    @staticmethod
    def _get_reldepth(gosubdag):
        """Get GO term depth from GoSubDag named tuples"""
        goid = next(iter(gosubdag.go_sources))
        return gosubdag.go2nt[goid].reldepth

    def chk_r0(self, arg):
        """Check relationship == False"""
        png, gosubdagplot = arg
        assert exists(png)
        go2p = gosubdagplot.gosubdag.rcntobj.go2ancestors
        exp_all = self.color2goids['green']
        ## print('AAAAAAAAAAAAAAAAAAAAAAAAA', exp_all)
        ## print('AAAAAAAAAAAAAAAAAAAAAAAAA', go2p[self.goid])
        self._cmp(go2p[self.goid], exp_all)
        assert not hasattr(self._get_source_term(gosubdagplot.gosubdag), 'reldepth')
        assert 'GO:0008150' not in go2p
        exp1 = {'GO:0008150'}  # biological process
        assert go2p['GO:0065007'] == exp1
        exp1.add('GO:0065007')
        self._cmp(go2p['GO:0050789'], exp1)
        exp1.add('GO:0050789')
        assert go2p['GO:0019222'] == exp1
        assert go2p['GO:0050794'] == exp1
        assert go2p['GO:0048518'] == exp1
        exp2 = set(exp1).union({'GO:0050794', 'GO:0048518'})
        assert go2p['GO:0048522'] == exp2
        #
        exp1.add('GO:0019222')
        assert go2p['GO:0060255'] == exp1
        exp1.add('GO:0060255')
        assert go2p['GO:0010468'] == exp1
        exp1.add('GO:0010468')
        exp1.add('GO:0050794')
        assert go2p['GO:0060968'] == exp1
        exp1.add('GO:0060968')
        assert go2p['GO:0060147'] == exp1
        exp1.add('GO:0060147')
        exp1.add('GO:0048522')
        assert go2p['GO:0060148'] == exp1.union(exp2)
        self._cmp(set(gosubdagplot.gosubdag.go2obj), exp_all.union({self.goid}))

    def _cmp(self, goset1, goset2):
        """Compare two sets of GO IDs"""
        assert goset1 == goset2, self._get_errstr(goset1, goset2)

    @staticmethod
    def _get_errstr(act_gos, exp_gos):
        """Get error string for mismatching sets of GO IDs"""
        return ' '.join([
            'ACT({A}) EXP({E})'.format(A=len(act_gos), E=len(exp_gos)),
            ' '.join(sorted(act_gos.symmetric_difference(exp_gos))),
        ])

    def _init_color2goids(self):
        """Initialize color-to-GO_IDs"""
        color2goids = defaultdict(set)
        for goid, hexcolor in self.go2color.items():
            color = self.hex2color[hexcolor]
            color2goids[color].add(goid)
        return color2goids

if __name__ == '__main__':
    test_getgo_color_for_plots()

# Copyright (C) 2016-present, DV Klopfenstein, H Tang, All rights reserved.
