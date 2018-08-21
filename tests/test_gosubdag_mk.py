#!/usr/bin/env python
"""Test all ways to make a GoSubDag."""

from goatools.gosubdag.gosubdag import GoSubDag
from goatools.test_data.nature3102_goea import get_goea_results

def test_gosubdag():
    """Test all ways to make a GoSubDag."""
    # Get data to use for test
    res = get_goea_results()  #
    godag = res['obo_dag']
    goea_results = res['goea_results']
    goids = [r.GO for r in goea_results]
    num_goids = len(goids)

    # Test Arg: goea_results (list of GOEnrichmentRecord objects)
    go_sources = [rec.GO for rec in goea_results]
    go2obj = {rec.GO:rec.goterm for rec in goea_results}
    gosubdag = GoSubDag(go_sources, go2obj)
    assert len(gosubdag.go_sources) == len(goea_results)

    # Test Arg: godag (GODag object)
    gosubdag = GoSubDag(None, godag)
    assert len(gosubdag.go_sources) > 40000

    # Test Arg: goids, godag
    gosubdag = GoSubDag(goids, godag)
    assert len(gosubdag.go_sources) == num_goids


if __name__ == '__main__':
    test_gosubdag()
