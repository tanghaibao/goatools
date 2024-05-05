#!/usr/bin/env python3
"""Small test utilities used by multiple tests"""

__copyright__ = (
    "Copyright (C) 2019-present, DV Klopfenstein, H Tang. All rights reserved."
)
__author__ = "DV Klopfenstein"

import sys
import timeit

from os.path import join
from os.path import dirname
from os.path import abspath

from datetime import timedelta
from goatools.base import get_godag as base_get_godag
from goatools.associations import dnld_annotation
from goatools.anno.factory import get_objanno as get_objanno_factory

from goatools.semantic import TermCounts


DIR_TEST = dirname(abspath(__file__))
REPO = abspath(join(DIR_TEST, ".."))


def prt_hms(tic, msg, prt=sys.stdout):
    """Print elapsed time and return current time"""
    toc = timeit.default_timer()
    prt.write("{HMS} {MSG}\n".format(HMS=str(timedelta(seconds=toc - tic)), MSG=msg))
    return toc


def repofn(fin):
    """Get a full filename, given a local file name from repo dir root"""
    return join(REPO, fin)


def get_godag(fin_godag, **kws):
    """Get GODAG containing only primary GO IDs (no alternate GO IDs)"""
    godag = base_get_godag(join(REPO, fin_godag), **kws)
    return {o.item_id: o for o in godag.values()}


def get_anno_fullname(fin_anno):
    """Get annotation filename"""
    fin_full = join(REPO, fin_anno)
    dnld_annotation(fin_full)
    return fin_full


def get_objanno(fin_anno, godag, namespace="all"):
    """Get annotation object"""
    fin_full = get_anno_fullname(fin_anno)
    return get_objanno_factory(fin_full, godag=godag, namespace=namespace)


def get_termcounts(fin_anno, godag, namespace="all", **kws):
    """Get termcounts object"""
    objanno = get_objanno(fin_anno, godag, namespace)
    id2gos = objanno.get_id2gos(namespace=namespace, **kws)
    return TermCounts(godag, id2gos)


# Copyright (C) 2019-present, DV Klopfenstein, et al. All rights reserved.
