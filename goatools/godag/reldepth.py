"""Manages a user-specified subset of a GO DAG."""

from __future__ import print_function

__copyright__ = "Copyright (C) 2016-present, DV Klopfenstein, H Tang, All rights reserved."
__author__ = "DV Klopfenstein"

from goatools.godag.consts import RELATIONSHIP_SET


def get_go2reldepth(goterms, relationship_set):
    """Calculate depth by tracing optional relationship_set"""
    return RelDepthCalc(relationship_set).get_go2reldepth_main(goterms)

# pylint: disable=too-few-public-methods
class RelDepthCalc:
    """Calculate depth by tracing optional relationship_set"""

    def __init__(self, relationship_set):
        self.relationship_set = relationship_set

    def get_go2reldepth_main(self, goterms):
        """Calculate depth by tracing optional relationship_set"""
        if self.relationship_set == RELATIONSHIP_SET:
            return self._get_go2reldepth_rall(goterms)
        return self._get_go2reldepth_rsome(goterms)

    @staticmethod
    def _get_go2reldepth_rall(goterms):
        """Calculate depth by tracing all optional relationship_set"""

        def _rall(rec, go2reldepth):
            """Get reldepth by tracing all optional relationship_set using recursive function"""
            if rec.item_id not in go2reldepth:
                anc = rec.get_goterms_upper()
                reldepth = max(_rall(rec, go2reldepth) for rec in anc) + 1 if anc else 0
                go2reldepth[rec.item_id] = reldepth
            return go2reldepth[rec.item_id]

        go2reldepth = {}
        for rec in goterms:
            if rec.item_id not in go2reldepth:
                _rall(rec, go2reldepth)
        return go2reldepth

    def _get_go2reldepth_rsome(self, goterms):
        """Calculate depth by tracing some optional relationship_set"""
        go2reldepth = {}
        for rec in goterms:
            if rec.item_id not in go2reldepth:
                self._rsome(rec, go2reldepth)
        return go2reldepth

    def _rsome(self, rec, go2reldepth):
        """Get reldepth by tracing some optional relationship_set using recursive function"""
        if rec.item_id not in go2reldepth:
            anc = rec.get_goterms_upper_rels(self.relationship_set)
            reldepth = max(self._rsome(rec, go2reldepth) for rec in anc) + 1 if anc else 0
            go2reldepth[rec.item_id] = reldepth
        return go2reldepth[rec.item_id]


# Copyright (C) 2016-present, DV Klopfenstein, H Tang, All rights reserved.
