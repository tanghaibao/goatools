"""Test fold_enrichment property on GOEnrichmentRecord."""

from goatools.go_enrichment import GOEnrichmentRecord

__copyright__ = "Copyright (C) 2010-2019, H Tang et al., All rights reserved."


def test_fold_enrichment_zero_pop_count():
    """fold_enrichment is None when pop_count is zero (term absent from population)."""
    rec = GOEnrichmentRecord(
        "GO:0000001",
        ratio_in_study=(5, 100),
        ratio_in_pop=(0, 1000),
        p_uncorrected=0.01,
        study_items=set(),
        pop_items=set(),
    )
    assert rec.fold_enrichment is None


def test_fold_enrichment_zero_study_n():
    """fold_enrichment is None when study_n is zero (empty study)."""
    rec = GOEnrichmentRecord(
        "GO:0000001",
        ratio_in_study=(0, 0),
        ratio_in_pop=(10, 1000),
        p_uncorrected=0.5,
        study_items=set(),
        pop_items=set(),
    )
    assert rec.fold_enrichment is None


def test_fold_enrichment_enriched():
    """fold_enrichment is > 1 for an enriched GO term and matches expected formula."""
    rec = GOEnrichmentRecord(
        "GO:0000001",
        ratio_in_study=(10, 100),
        ratio_in_pop=(10, 1000),
        p_uncorrected=0.001,
        study_items=set(),
        pop_items=set(),
    )
    fe = rec.fold_enrichment
    # (10/100) / (10/1000) = 0.1 / 0.01 = 10.0
    assert fe is not None
    assert abs(fe - 10.0) < 1e-10
    assert rec.enrichment == "e"
    assert fe > 1.0


def test_fold_enrichment_purified():
    """fold_enrichment is < 1 for a purified GO term and matches expected formula."""
    rec = GOEnrichmentRecord(
        "GO:0000001",
        ratio_in_study=(1, 100),
        ratio_in_pop=(100, 1000),
        p_uncorrected=0.001,
        study_items=set(),
        pop_items=set(),
    )
    fe = rec.fold_enrichment
    # (1/100) / (100/1000) = 0.01 / 0.1 = 0.1
    assert fe is not None
    assert abs(fe - 0.1) < 1e-10
    assert rec.enrichment == "p"
    assert fe < 1.0


def test_fold_enrichment_equal():
    """fold_enrichment is 1.0 when study and population proportions are equal."""
    rec = GOEnrichmentRecord(
        "GO:0000001",
        ratio_in_study=(10, 100),
        ratio_in_pop=(100, 1000),
        p_uncorrected=1.0,
        study_items=set(),
        pop_items=set(),
    )
    fe = rec.fold_enrichment
    # (10/100) / (100/1000) = 0.1 / 0.1 = 1.0
    assert fe is not None
    assert abs(fe - 1.0) < 1e-10


def test_fold_enrichment_zero_study_count():
    """fold_enrichment is 0.0 when study_count is zero but study_n is non-zero."""
    rec = GOEnrichmentRecord(
        "GO:0000001",
        ratio_in_study=(0, 100),
        ratio_in_pop=(50, 1000),
        p_uncorrected=0.5,
        study_items=set(),
        pop_items=set(),
    )
    fe = rec.fold_enrichment
    # (0/100) / (50/1000) = 0.0 / 0.05 = 0.0
    assert fe is not None
    assert abs(fe - 0.0) < 1e-10
    assert rec.enrichment == "p"


if __name__ == "__main__":
    test_fold_enrichment_zero_pop_count()
    test_fold_enrichment_zero_study_n()
    test_fold_enrichment_enriched()
    test_fold_enrichment_purified()
    test_fold_enrichment_equal()
    print("All fold_enrichment tests passed.")
