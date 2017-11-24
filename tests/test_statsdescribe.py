#!/usr/bin/env python
"""Test StatsDescribe."""

from goatools.statsdescribe import StatsDescribe
from goatools.test_data.gjoneska_goea_consistent_increase import goea_results

def test_statsdescribe():
    """Use StatsDescribe to create a markdown table.

fdr_bh
name     | # fdr_bh | range of fdr_bh      | 25th perc|   median | 75th perc|     mean | stddev
---------|----------|----------------------|----------|----------|----------|----------|---------
GOATOOLS |       59 | 1.87e-07 to 4.94e-02 | 2.72e-04 | 1.03e-02 | 3.04e-02 | 1.56e-02 | 1.82e-02

    """
    #pylint: disable=no-member
    # Somehow goea_results contains fields of empty string, which we can check with:
    # print([(nt.GO, nt.p_fdr_bh) for nt in goea_results])
    nts_goids = [nt for nt in goea_results if nt.p_fdr_bh != '' and nt.p_fdr_bh < 0.05]
    fdr_vals = [nt.p_fdr_bh for nt in nts_goids]
    statsobj = StatsDescribe("fdr_bh", fmtstr="{:>8.2e}")
    statsobj.prt_hdr()
    statsobj.prt_data("GOATOOLS", fdr_vals)


if __name__ == '__main__':
    test_statsdescribe()
