#!/usr/bin/env python
"""Tests that GafReader supports use_symbol option to key annotations by DB_Symbol.

This tests the fix for the issue where users with gene symbols in their population
and study files get "FATAL: NO POPULATION ITEMS SEEN IN THE ANNOTATIONS" when using
a GAF file that stores UniProt IDs in DB_ID and gene symbols in DB_Symbol.
"""

import os
from goatools.anno.gaf_reader import GafReader
from goatools.base import get_godag

REPO = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..")
FIN_GAF = os.path.join(REPO, "tests/data/gaf_with_uniprot_ids.gaf")
FIN_OBO = os.path.join(REPO, "tests/data/goslim_generic.obo")


def test_gaf_use_symbol_default():
    """Default behavior: annotations keyed by DB_ID (UniProt accessions)."""
    gaf = GafReader(FIN_GAF, prt=None)
    id2gos = gaf.get_id2gos(namespace="BP", prt=None)
    # Keys should be DB_ID values (UniProt IDs like P12345)
    assert id2gos, "Expected non-empty annotations"
    for key in id2gos:
        assert key.startswith("P"), f"Expected UniProt ID key, got: {key}"
    # Gene symbols should NOT be present as keys
    assert "GENE1" not in id2gos, "Gene symbol GENE1 should not be a key in default mode"
    print(f"Default mode: {len(id2gos)} IDs (UniProt) in BP annotations")


def test_gaf_use_symbol_true():
    """use_symbol=True: annotations keyed by DB_Symbol (gene symbols)."""
    gaf = GafReader(FIN_GAF, prt=None, use_symbol=True)
    id2gos = gaf.get_id2gos(namespace="BP", prt=None)
    # Keys should be DB_Symbol values (gene symbols like GENE1)
    assert id2gos, "Expected non-empty annotations"
    for key in id2gos:
        assert key.startswith("GENE"), f"Expected gene symbol key, got: {key}"
    # UniProt IDs should NOT be present as keys
    assert "P12345" not in id2gos, "UniProt ID P12345 should not be a key in use_symbol mode"
    print(f"use_symbol mode: {len(id2gos)} gene symbols in BP annotations")


def test_gaf_use_symbol_propagate_counts():
    """use_symbol=True with propagate_counts: annotations keyed by DB_Symbol."""
    godag = get_godag(FIN_OBO, optional_attrs=set())
    gaf = GafReader(FIN_GAF, prt=None, use_symbol=True, godag=godag)
    id2gos = gaf.get_id2gos(namespace="BP", prt=None)
    # First, verify use_symbol=True with propagate_counts works
    # (calls _get_dbid2goids_p1)
    id2gos_prop = gaf._get_id2gos(
        [nt for nt in gaf.associations if nt.NS == "BP"],
        propagate_counts=True,
        prt=None,
    )
    assert id2gos_prop, "Expected non-empty annotations with propagate_counts"
    for key in id2gos_prop:
        assert key.startswith("GENE"), f"Expected gene symbol key with propagate, got: {key}"
    print(f"use_symbol+propagate mode: {len(id2gos_prop)} gene symbols in BP annotations")


def test_gaf_use_symbol_goea():
    """Test that GOEA works when using gene symbols as population/study IDs."""
    from goatools.goea.go_enrichment_ns import GOEnrichmentStudyNS

    godag = get_godag(FIN_OBO, optional_attrs=set())
    gaf = GafReader(FIN_GAF, prt=None, use_symbol=True, godag=godag)

    ns2assc = gaf.get_ns2assc()
    # Population: all gene symbols in annotations
    pop = set()
    for ids in ns2assc.values():
        pop.update(ids.keys())
    assert pop, "Population should be non-empty"
    # Study: a small subset of gene symbols
    study = frozenset(["GENE1", "GENE2", "GENE3"])
    assert study.issubset(pop), "Study genes should be in the population"

    objgoeans = GOEnrichmentStudyNS(
        pop,
        ns2assc,
        godag,
        propagate_counts=False,
        alpha=0.05,
        methods=["bonferroni"],
    )
    results = objgoeans.run_study(study)
    assert results is not None, "GOEA results should not be None"
    print(f"GOEA with use_symbol: {len(results)} results")
    print("TEST PASSED")


if __name__ == "__main__":
    test_gaf_use_symbol_default()
    test_gaf_use_symbol_true()
    test_gaf_use_symbol_propagate_counts()
    test_gaf_use_symbol_goea()
