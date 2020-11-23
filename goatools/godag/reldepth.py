"""Manages a user-specified subset of a GO DAG."""

from __future__ import print_function

__copyright__ = "Copyright (C) 2016-present, DV Klopfenstein, H Tang, All rights reserved."
__author__ = "DV Klopfenstein"

from goatools.godag.consts import RELATIONSHIP_SET


def get_go2reldepth(godag, relationships):
    """Calculate depth by tracing optional relationships"""
    return RelDepthCalc(godag, relationships).get_go2reldepth()

# pylint: disable=too-few-public-methods
class RelDepthCalc:
    """Calculate depth by tracing optional relationships"""

    def __init__(self, godag, relationships):
        self.godag = godag
        self.relationships = relationships

    def get_go2reldepth(self):
        """Calculate depth by tracing optional relationships"""
        if self.relationships == RELATIONSHIP_SET:
            return self._get_go2reldepth_rall()
        return self._get_go2reldepth_rsome()

    def _get_go2reldepth_rall(self):
        """Calculate depth by tracing all optional relationships"""
        go2reldepth = {}
        for rec in self.godag.values():
            self._rall(rec, go2reldepth)
        return go2reldepth

    def _get_go2reldepth_rsome(self):
        """Calculate depth by tracing some optional relationships"""
        go2reldepth = {}
        for rec in self.godag.values():
            self._rsome(rec, go2reldepth)
        return go2reldepth

    def _rall(self, rec, go2reldepth):
        """Get reldepth by tracing all optional relationships using recursive function"""
        if rec.item_id not in go2reldepth:
            anc = rec.get_goterms_upper()
            reldepth = max(self._rall(rec, go2reldepth) for rec in anc) + 1 if anc else 0
            go2reldepth[rec.item_id] = reldepth
        return go2reldepth[rec.item_id]

    def _rsome(self, rec, go2reldepth):
        """Get reldepth by tracing some optional relationships using recursive function"""
        if rec.item_id not in go2reldepth:
            anc = rec.get_goterms_upper_rels(self.relationships)
            reldepth = max(self._rsome(rec, go2reldepth) for rec in anc) + 1 if anc else 0
            go2reldepth[rec.item_id] = reldepth
        return go2reldepth[rec.item_id]


# Copyright (C) 2016-present, DV Klopfenstein, H Tang, All rights reserved.
