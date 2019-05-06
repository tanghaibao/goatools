"""Keyword args for selecting annotation lines."""

__copyright__ = "Copyright (C) 2016-2019, DV Klopfenstein, H Tang. All rights reserved."
__author__ = "DV Klopfenstein"


# pylint: disable=too-few-public-methods
class AnnoOptions(object):
    """Keyword args for selecting annotation lines."""

    keys_exp = set(['ev_include', 'ev_exclude',
                    'b_geneid2gos', 'go2geneids',
                    'namespace', 'keep_ND', 'keep_NOT'])

    def __init__(self, evobj, **kws):
        # Get associations only for specified Evidence_Codes
        # kws: ev_include ev_exclude keep_ND keep_NOT b_geneid2gos go2geneids
        self.evobj = evobj
        _incexc2codes = evobj.get_min_inc_exc(kws.get('ev_include'), kws.get('ev_exclude'))
        self.include_evcodes = _incexc2codes.get('inc')
        self.exclude_evcodes = _incexc2codes.get('exc')
        # Associations are normally gene2gos
        # Return go2genes, if desired
        self.b_geneid2gos = self._init_b_geneid2gos(kws)
        # Keep annotation, even if:
        #   * Evidence_Code == ND -> No biological data No biological Data available
        self._keep_nd = kws.get('keep_ND', False)
        #   * Qualifiers contain NOT
        self._keep_not = kws.get('keep_NOT', False)


        # keep_qualified keep_unqualified keep_nd keep_not ev_include ev_exclude
        self.param2fnc = {
            ('keep_qualified', 0): self._qual_0,
            ('keep_qualified', 1): self._qual_inc,
            ('keep_qualified', -1): self._qual_exc,

            ('keep_unqualified', 0): self._all_0,
            ('keep_unqualified', 1): self._all_inc,
            ('keep_unqualified', -1): self._all_exc,

            ('keep_ND', 0): self._nd_0,
            ('keep_ND', 1): self._nd_inc,
            ('keep_ND', -1): self._nd_exc,

            ('keep_NOT', 0): self._not_0,
            ('keep_NOT', 1): self._not_inc,
            ('keep_NOT', -1): self._not_exc,
        }

    # pylint: disable=bad-whitespace
    nd_not2desc = {
        # Keep ND  Keep NOT
        (False, False): 'keep_qualified',
        (True,  True):  'keep_unqualified',
        (True,  False): 'keep_ND',
        (False, True):  'keep_NOT',
    }

    incexc2num = {
        (False, False): 0,  # Default evidence codes
        (True,  False): 1,  # Include user-specified Evidence codes
        (False, True): -1,  # Exclude user-specified Evidence codes
    }

    def __str__(self):
        pat = (
            'EVIDENCE CODE ARGS: '
            '{DESC}(ND={ND:1}, NOT={NOT:1}) CODES(inc[{I}], exc[{E}]) '
            '{A2Bs}'
        )
        return pat.format(
            DESC=self._get_desc(), ND=self._keep_nd, NOT=self._keep_not,
            I='-' if self.include_evcodes is None else len(self.include_evcodes),
            E='-' if self.exclude_evcodes is None else len(self.exclude_evcodes),
            A2Bs='id2gos' if self.b_geneid2gos else 'go2ids',
        )

    def getfnc_qual_ev(self):
        """Keep annotaion if it passes potentially modified selection."""
        fnc_key = (
            self._get_desc(),
            self.incexc2num[(
                self.include_evcodes is not None,
                self.exclude_evcodes is not None)],
        )
        return self.param2fnc[fnc_key]

    def _get_desc(self):
        """Get desc: qualified, unqualified; ND and NOT"""
        return self.nd_not2desc[(self._keep_nd, self._keep_not)]


    # - Filter by Evidence_code or Qualifier ----------------------------------------
    # pylint: disable=unused-argument
    @staticmethod
    def _qual_0(qualifiers, evidence_code):
        """Keep qualified; NO evidence filter"""
        return 'NOT' not in qualifiers and evidence_code != 'ND'

    def _qual_inc(self, qualifiers, evidence_code):
        """Keep qualified; IN evidence filter"""
        return 'NOT' not in qualifiers and evidence_code != 'ND' and \
            evidence_code in self.include_evcodes

    def _qual_exc(self, qualifiers, evidence_code):
        """Keep qualified; IN evidence filter"""
        return 'NOT' not in qualifiers and evidence_code != 'ND' and \
            evidence_code not in self.exclude_evcodes

    @staticmethod
    def _all_0(qualifiers, evidence_code):
        """Keep unqualified: NOT in Qualifiers AND ND in Evidence code; NO evidence filter"""
        return True

    def _all_inc(self, qualifiers, evidence_code):
        """Keep unqualified: NOT in Qualifier AND ND in Evidence code; IN evidence filter"""
        return evidence_code in self.include_evcodes

    def _all_exc(self, qualifiers, evidence_code):
        """Keep unqualified: NOT in Qualifier AND ND in Evidence code; IN evidence filter"""
        return evidence_code not in self.exclude_evcodes

    @staticmethod
    def _nd_0(qualifiers, evidence_code):
        """Keep ND; NO evidence filter"""
        return 'NOT' not in qualifiers

    def _nd_inc(self, qualifiers, evidence_code):
        """Keep ND; IN evidence filter"""
        return 'NOT' not in qualifiers and evidence_code in self.include_evcodes

    def _nd_exc(self, qualifiers, evidence_code):
        """Keep ND; IN evidence filter"""
        return 'NOT' not in qualifiers and evidence_code not in self.exclude_evcodes

    @staticmethod
    def _not_0(qualifiers, evidence_code):
        """Keep NOT; NO evidence filter"""
        return evidence_code != 'ND'

    def _not_inc(self, qualifiers, evidence_code):
        """Keep NOT; IN evidence filter"""
        return evidence_code != 'ND' and evidence_code in self.include_evcodes

    def _not_exc(self, qualifiers, evidence_code):
        """Keep NOT; IN evidence filter"""
        return evidence_code != 'ND' and evidence_code not in self.exclude_evcodes

    # -------------------------------------------------------------------------------
    @staticmethod
    def _init_b_geneid2gos(kws):
        """GEt b_geneid2gos."""
        if 'b_geneid2gos' in kws:
            return kws['b_geneid2gos']
        if 'go2geneids' in kws:
            return not kws['go2geneids']
        return True


# Copyright (C) 2016-2019, DV Klopfenstein, H Tang. All rights reserved."
