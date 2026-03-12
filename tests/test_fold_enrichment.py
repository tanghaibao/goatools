"""Test fold_enrichment property on GOEnrichmentRecord."""

from goatools.go_enrichment import GOEnrichmentRecord

__copyright__ = "Copyright (C) 2010-2019, H Tang et al., All rights reserved."


def _make_rec(ratio_in_study, ratio_in_pop, p_uncorrected=0.01, depth=5):
    """Helper: create a GOEnrichmentRecord with given ratios."""
    return GOEnrichmentRecord(
        "GO:0000001",
        ratio_in_study=ratio_in_study,
        ratio_in_pop=ratio_in_pop,
        p_uncorrected=p_uncorrected,
        study_items=set(),
        pop_items=set(),
        depth=depth,
    )


def test_fold_enrichment_zero_pop_count():
    """fold_enrichment is None when pop_count is zero (term absent from population)."""
    rec = _make_rec(ratio_in_study=(5, 100), ratio_in_pop=(0, 1000))
    assert rec.fold_enrichment is None


def test_fold_enrichment_zero_study_n():
    """fold_enrichment is None when study_n is zero (empty study)."""
    rec = _make_rec(ratio_in_study=(0, 0), ratio_in_pop=(10, 1000), p_uncorrected=0.5)
    assert rec.fold_enrichment is None


def test_fold_enrichment_enriched():
    """fold_enrichment is > 1 for an enriched GO term and matches expected formula."""
    rec = _make_rec(ratio_in_study=(10, 100), ratio_in_pop=(10, 1000), p_uncorrected=0.001)
    fe = rec.fold_enrichment
    # (10/100) / (10/1000) = 0.1 / 0.01 = 10.0
    assert fe is not None
    assert abs(fe - 10.0) < 1e-10
    assert rec.enrichment == "e"
    assert fe > 1.0


def test_fold_enrichment_purified():
    """fold_enrichment is < 1 for a purified GO term and matches expected formula."""
    rec = _make_rec(ratio_in_study=(1, 100), ratio_in_pop=(100, 1000), p_uncorrected=0.001)
    fe = rec.fold_enrichment
    # (1/100) / (100/1000) = 0.01 / 0.1 = 0.1
    assert fe is not None
    assert abs(fe - 0.1) < 1e-10
    assert rec.enrichment == "p"
    assert fe < 1.0


def test_fold_enrichment_equal():
    """fold_enrichment is 1.0 when study and population proportions are equal."""
    rec = _make_rec(ratio_in_study=(10, 100), ratio_in_pop=(100, 1000), p_uncorrected=1.0)
    fe = rec.fold_enrichment
    # (10/100) / (100/1000) = 0.1 / 0.1 = 1.0
    assert fe is not None
    assert abs(fe - 1.0) < 1e-10


def test_fold_enrichment_zero_study_count():
    """fold_enrichment is 0.0 when study_count is zero but study_n is non-zero."""
    rec = _make_rec(ratio_in_study=(0, 100), ratio_in_pop=(50, 1000), p_uncorrected=0.5)
    fe = rec.fold_enrichment
    # (0/100) / (50/1000) = 0.0 / 0.05 = 0.0
    assert fe is not None
    assert abs(fe - 0.0) < 1e-10
    assert rec.enrichment == "p"


def test_fold_enrichment_in_default_fields():
    """fold_enrichment is included in the default print fields."""
    assert "fold_enrichment" in GOEnrichmentRecord._fldsdefprt


def test_fold_enrichment_in_prtflds_default():
    """fold_enrichment is included in get_prtflds_default() output."""
    rec = _make_rec(ratio_in_study=(10, 100), ratio_in_pop=(10, 1000))
    assert "fold_enrichment" in rec.get_prtflds_default()


def test_str_with_valid_fold_enrichment():
    """__str__ includes fold_enrichment value when it is a valid float."""
    rec = _make_rec(ratio_in_study=(10, 100), ratio_in_pop=(10, 1000), p_uncorrected=0.001)
    result = str(rec)
    fields = result.split("\t")
    # fold_enrichment is at index 6 (7th field) in _fldsdefprt
    fe_idx = GOEnrichmentRecord._fldsdefprt.index("fold_enrichment")
    assert abs(float(fields[fe_idx]) - 10.0) < 1e-3


def test_str_with_none_fold_enrichment():
    """__str__ shows 'n.a.' for fold_enrichment when it cannot be computed (pop_count=0)."""
    rec = _make_rec(ratio_in_study=(5, 100), ratio_in_pop=(0, 1000))
    result = str(rec)
    fields = result.split("\t")
    fe_idx = GOEnrichmentRecord._fldsdefprt.index("fold_enrichment")
    assert fields[fe_idx] == "n.a."


def test_get_field_values_with_none_fold_enrichment():
    """get_field_values returns 'n.a.' for fold_enrichment when it is None."""
    rec = _make_rec(ratio_in_study=(5, 100), ratio_in_pop=(0, 1000))
    row = rec.get_field_values(["fold_enrichment"], rpt_fmt=True)
    assert row == ["n.a."]


def test_get_field_values_fold_enrichment_raw():
    """get_field_values returns the raw None when rpt_fmt=False and fold_enrichment is None."""
    rec = _make_rec(ratio_in_study=(5, 100), ratio_in_pop=(0, 1000))
    row = rec.get_field_values(["fold_enrichment"], rpt_fmt=False)
    assert row == [None]


def test_get_field_values_fold_enrichment_float():
    """get_field_values returns the float fold_enrichment value when it is valid."""
    rec = _make_rec(ratio_in_study=(10, 100), ratio_in_pop=(10, 1000))
    row = rec.get_field_values(["fold_enrichment"], rpt_fmt=False)
    assert len(row) == 1
    assert abs(row[0] - 10.0) < 1e-10


if __name__ == "__main__":
    test_fold_enrichment_zero_pop_count()
    test_fold_enrichment_zero_study_n()
    test_fold_enrichment_enriched()
    test_fold_enrichment_purified()
    test_fold_enrichment_equal()
    test_fold_enrichment_in_default_fields()
    test_fold_enrichment_in_prtflds_default()
    test_str_with_valid_fold_enrichment()
    test_str_with_none_fold_enrichment()
    test_get_field_values_with_none_fold_enrichment()
    test_get_field_values_fold_enrichment_raw()
    test_get_field_values_fold_enrichment_float()
    print("All fold_enrichment tests passed.")
