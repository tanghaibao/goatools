"""Keyword args for selecting annotation lines."""

__copyright__ = "Copyright (C) 2016-2019, DV Klopfenstein, H Tang. All rights reserved."
__author__ = "DV Klopfenstein"


class AnnoOptions:
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

    def keep(self, qualifiers, evidence_code):
        """Keep annotaion if it passes potentially modified selection."""
        return self.keep_qualified(qualifiers, evidence_code) and self.keep_evidence(evidence_code)

    def keep_qualified(self, qualifiers, evidence_code):
        """Normally keeps qualified associations, but can keep more."""
        # NOT: Used when gene is expected to have function F, but does NOT.
        # ND : GO function not seen after exhaustive annotation attempts to the gene.
        if self._keep_qualified:
            # Keep everything but these:
            #     Qualifiers contain NOT
            #     Evidence_Code == ND -> No biological data No biological Data available
            return 'not' not in qualifiers and evidence_code != 'ND'
        if self._keep_unqualified:
            return True
        if self._keep_nd:
            return 'not' not in qualifiers
        if self._keep_not:
            return evidence_code != 'ND'

    def keep_evidence(self, evidence_code):
        """Keep all evidence (default) or selected Evidence_Codes."""
        if self.evidence_set is None:
            return True
        return evidence_code in self.evidence_set

    @staticmethod
    def _init_b_geneid2gos(kws):
        """GEt b_geneid2gos."""
        if 'b_geneid2gos' in kws:
            return kws['b_geneid2gos']
        if 'go2geneids' in kws:
            return not kws['go2geneids']
        return True


# Copyright (C) 2016-2019, DV Klopfenstein, H Tang. All rights reserved."
