#!/usr/bin/env python
"""_resultsTest writing GOATOOLS Gene Ontology Gene Enrichment results to a Python module."""

from __future__ import print_function

from goatools.test_data.nature3102_goea import get_goea_results
from goatools.rpt.goea_nt_xfrm import get_goea_nts_prt
from goatools.rpt.goea_nt_xfrm import MgrNtGOEAs


def test_wrpy():
    """Test writing GOATOOLS GOEA results to a Python module as a list of nts."""
    # Run GOATOOLS Gene Ontology Enrichment Analysis
    nature_data = get_goea_results()
    # Convert GOATOOLS GOEA results into a list of namedtuples
    goea_results = [nt for nt in nature_data['goea_results'] if nt.p_fdr_bh < 0.05]
    nts_goea = get_goea_nts_prt(goea_results)
    _run(nts_goea, next(iter(nts_goea))._fields)
    # print(vars(next(iter(goea_results))).keys())

def _run(goea_results, flds):
    """Get string representations of P-value data floats ."""
    obj = MgrNtGOEAs(goea_results)
    pval_flds = set("s_" + k[2:] for k in flds if k[:2] == "p_")
    nts = obj.get_nts_strpval()
    fmts = ["{FLD}={{{FLD}}}".format(FLD=f) for f in pval_flds]
    fmt = "{N} {NS} {GO} " + " ".join(fmts) + " {name}"
    for idx, ntd in enumerate(nts):
        print(fmt.format(N=idx, **ntd._asdict()))
    print("FORMAT", fmt)


if __name__ == '__main__':
    test_wrpy()
