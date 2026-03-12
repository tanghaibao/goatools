"""Test filtering enrichment analysis to user-specified GO terms."""

# pylint: disable=invalid-name

import os

from goatools.associations import read_associations
from goatools.go_enrichment import GOEnrichmentStudy
from goatools.obo_parser import GODag

__copyright__ = "Copyright (C) 2010-2024, DV Klopfenstein, H Tang. All rights reserved."

REPO = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..")


def _get_goeaobj(methods=None):
    """Create a GOEnrichmentStudy with local test data."""
    obo_fin = os.path.join(REPO, "data/i86.obo")
    obo_dag = GODag(obo_fin, load_obsolete=True)
    fin_assc = os.path.join(REPO, "tests/data/small_association")
    assoc = read_associations(fin_assc, "id2gos", no_top=True)
    popul_fin = os.path.join(REPO, "tests/data/small_population")
    popul_ids = [line.rstrip() for line in open(popul_fin, encoding="utf-8")]
    goeaobj = GOEnrichmentStudy(popul_ids, assoc, obo_dag, methods=methods)
    return goeaobj


def test_goea_selected_goids():
    """Test that selected_goids filters results to only the specified GO terms."""
    goeaobj = _get_goeaobj(methods=["bonferroni"])
    study_fin = os.path.join(REPO, "tests/data/small_study")
    study_ids = [line.rstrip() for line in open(study_fin, encoding="utf-8")]

    # Run without filter to get baseline
    results_all = goeaobj.run_study(study_ids, log=None)
    go_ids_all = set(r.GO for r in results_all)
    assert len(go_ids_all) > 2, "Expected more than 2 GO terms in unfiltered results"

    # Run with a subset of GO terms
    selected = {"GO:0008150", "GO:0003674"}
    results_filtered = goeaobj.run_study(study_ids, selected_goids=selected, log=None)
    go_ids_filtered = set(r.GO for r in results_filtered)

    # Verify only the requested GO terms are in the results
    assert go_ids_filtered.issubset(selected), (
        "Filtered results contain GO terms not in selected_goids: "
        "{EXTRA}".format(EXTRA=go_ids_filtered - selected)
    )
    assert len(results_filtered) < len(results_all), (
        "Filtered results should be fewer than unfiltered results"
    )
    print("PASSED test_goea_selected_goids: "
          "{N} filtered results (from {M} total)".format(
              N=len(results_filtered), M=len(results_all)))


def test_goea_selected_goids_single():
    """Test filtering to a single GO term."""
    goeaobj = _get_goeaobj(methods=["bonferroni"])
    study_fin = os.path.join(REPO, "tests/data/small_study")
    study_ids = [line.rstrip() for line in open(study_fin, encoding="utf-8")]

    selected = {"GO:0008150"}
    results_filtered = goeaobj.run_study(study_ids, selected_goids=selected, log=None)
    go_ids_filtered = set(r.GO for r in results_filtered)

    assert go_ids_filtered.issubset(selected), (
        "Filtered results contain unexpected GO terms: "
        "{EXTRA}".format(EXTRA=go_ids_filtered - selected)
    )
    print("PASSED test_goea_selected_goids_single: "
          "{N} result(s) for 1 selected GO term".format(N=len(results_filtered)))


def test_goea_selected_goids_none():
    """Test that passing selected_goids=None returns all GO terms (default behavior)."""
    goeaobj = _get_goeaobj(methods=["bonferroni"])
    study_fin = os.path.join(REPO, "tests/data/small_study")
    study_ids = [line.rstrip() for line in open(study_fin, encoding="utf-8")]

    results_all = goeaobj.run_study(study_ids, log=None)
    results_none_filter = goeaobj.run_study(study_ids, selected_goids=None, log=None)

    assert len(results_all) == len(results_none_filter), (
        "selected_goids=None should return the same number of results as no filter"
    )
    print("PASSED test_goea_selected_goids_none: "
          "{N} results with and without filter".format(N=len(results_all)))


if __name__ == "__main__":
    test_goea_selected_goids()
    test_goea_selected_goids_single()
    test_goea_selected_goids_none()

# Copyright (C) 2010-2024, DV Klopfenstein, H Tang. All rights reserved.
