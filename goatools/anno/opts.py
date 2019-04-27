"""Keyword args for selecting annotation lines."""

__copyright__ = "Copyright (C) 2016-2019, DV Klopfenstein, H Tang. All rights reserved."
__author__ = "DV Klopfenstein"


class AnnoOptions(object):
    """Keyword args for selecting annotation lines."""

    keys_exp = set(['evidence_set',
                    'b_geneid2gos', 'go2geneids',
                    'keep_ND', 'keep_NOT'])

    def __init__(self, **kws):
        # Get associations only for specified Evidence_Codes
        # kws: evidence_set keep_ND keep_NOT b_geneid2gos go2geneids
        self.evidence_set = kws.get('evidence_set', None)
        # Associations are normally gene2gos
        # Return go2genes, if desired
        self.b_geneid2gos = self._init_b_geneid2gos(kws)
        # Keep annotation, even if:
        #   * Evidence_Code == ND -> No biological data No biological Data available
        self._keep_nd = kws.get('keep_ND', False)
        #   * Qualifiers contain NOT
        self._keep_not = kws.get('keep_NOT', False)
        self._keep_qualified = not self._keep_nd and not self._keep_not
        self._keep_unqualified = self._keep_nd and self._keep_not

        # keep_qualified keep_unqualified keep_nd keep_not evidence_set
        # pylint: disable=bad-whitespace
        self.param2fnc = {
            ('keep_qualified',    True):  self._qual_0,
            ('keep_qualified',    False): self._qual_1,
            ('keep_unqualified',  True):  self._all_0,
            ('keep_unqualified',  False): self._all_1,
            ('keep_nd',           True):  self._nd_0,
            ('keep_nd',           False): self._nd_1,
            ('keep_not',          True):  self._not_0,
            ('keep_not',          False): self._not_0,
        }

    def getfnc_qual_ev(self):
        """Keep annotaion if it passes potentially modified selection."""
        fnc_key = (
            self._get_keep_key(),
            self.evidence_set is None,
        )
        return self.param2fnc[fnc_key]

    def _get_keep_key(self):
        """Normally keeps qualified associations, but can keep more."""
        # NOT: Used when gene is expected to have function F, but does NOT.
        # ND : GO function not seen after exhaustive annotation attempts to the gene.
        # print('------------------------- AnnoOptions::keep_qualified(', qualifiers, evidence_code)
        if self._keep_qualified:
            # Keep everything but these:
            #     Qualifiers contain NOT
            #     Evidence_Code == ND -> No biological data No biological Data available
            return 'keep_qualified'
        if self._keep_unqualified:
            return 'keep_unqualified'
        if self._keep_nd:
            return 'keep_nd'
        if self._keep_not:
            return 'keep_not'

    # - Filter by Evidence_code or Qualifier ----------------------------------------
    # pylint: disable=unused-argument
    @staticmethod
    def _qual_0(qualifiers, evidence_code):
        """Keep qualified; NO evidence filter"""
        return 'NOT' not in qualifiers and evidence_code != 'ND'

    def _qual_1(self, qualifiers, evidence_code):
        """Keep qualified; IN evidence filter"""
        return 'NOT' not in qualifiers and evidence_code != 'ND' and \
            evidence_code in self.evidence_set

    @staticmethod
    def _all_0(qualifiers, evidence_code):
        """Keep unqualified: NOT in Qualifiers AND ND in Evidence code; NO evidence filter"""
        return True

    def _all_1(self, qualifiers, evidence_code):
        """Keep unqualified: NOT in Qualifier AND ND in Evidence code; IN evidence filter"""
        return evidence_code in self.evidence_set

    @staticmethod
    def _nd_0(qualifiers, evidence_code):
        """Keep ND; NO evidence filter"""
        return 'NOT' not in qualifiers

    def _nd_1(self, qualifiers, evidence_code):
        """Keep ND; IN evidence filter"""
        return 'NOT' not in qualifiers and evidence_code in self.evidence_set

    @staticmethod
    def _not_0(qualifiers, evidence_code):
        """Keep NOT; NO evidence filter"""
        return evidence_code != 'ND'

    def _not_1(self, qualifiers, evidence_code):
        """Keep NOT; IN evidence filter"""
        return evidence_code != 'ND' and evidence_code in self.evidence_set

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
