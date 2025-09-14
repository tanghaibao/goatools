#!/usr/bin/env python
"""Comprehensive test cases for goatools.semantic module.

This module provides extensive test coverage for all functions and classes
in the semantic.py module, including edge cases and error conditions.
"""

import pytest

from unittest.mock import patch


from goatools.semantic import (
    TermCounts,
    get_info_content,
    resnik_sim,
    lin_sim,
    lin_sim_calc,
    schlicker_sim,
    schlicker_sim_calc,
    get_freq_msca,
    common_parent_go_ids,
    deepest_common_ancestor,
    min_branch_length,
    semantic_distance,
    semantic_similarity,
)
from tests.utils import get_godag


class TestTermCounts:
    """Test cases for the TermCounts class."""

    def setup_method(self):
        """Set up test fixtures."""
        self.godag = get_godag("go-basic.obo")
        self.annots = {
            "gene1": {"GO:0008150", "GO:0003674"},  # BP and MF
            "gene2": {"GO:0008150", "GO:0005575"},  # BP and CC
            "gene3": {"GO:0003674", "GO:0005575"},  # MF and CC
        }

    def test_init_basic(self):
        """Test basic TermCounts initialization."""
        termcounts = TermCounts(self.godag, self.annots)

        assert termcounts.go2obj == self.godag
        assert len(termcounts.gocnts) > 0
        assert len(termcounts.go2genes) > 0
        assert len(termcounts.gene2gos) == 3
        assert "GO:0008150" in termcounts.gocnts
        assert "GO:0003674" in termcounts.gocnts
        assert "GO:0005575" in termcounts.gocnts

    def test_init_with_relationships(self):
        """Test TermCounts initialization with specific relationships."""
        relationships = {"is_a", "part_of"}
        termcounts = TermCounts(self.godag, self.annots, relationships=relationships)

        assert termcounts.gosubdag is not None
        # When GODag is not loaded with relationships, an empty set is expected
        # This tests the correct behavior when relationships are requested but not available
        assert termcounts.gosubdag.relationships == set()

    def test_init_with_prt(self):
        """Test TermCounts initialization with print output."""
        with patch("sys.stdout"):
            termcounts = TermCounts(self.godag, self.annots, prt=True)
            # Should not raise an exception and should have printed something
            assert termcounts is not None

    def test_get_count(self):
        """Test getting count for a specific GO term."""
        termcounts = TermCounts(self.godag, self.annots)

        # Test existing GO term
        count = termcounts.get_count("GO:0008150")
        assert isinstance(count, int)
        assert count >= 0

        # Test non-existing GO term
        count = termcounts.get_count("GO:9999999")
        assert count == 0

    def test_get_total_count(self):
        """Test getting total count for each aspect."""
        termcounts = TermCounts(self.godag, self.annots)

        bp_count = termcounts.get_total_count("biological_process")
        mf_count = termcounts.get_total_count("molecular_function")
        cc_count = termcounts.get_total_count("cellular_component")

        assert isinstance(bp_count, int)
        assert isinstance(mf_count, int)
        assert isinstance(cc_count, int)
        assert bp_count >= 0
        assert mf_count >= 0
        assert cc_count >= 0

    def test_get_term_freq(self):
        """Test getting term frequency."""
        termcounts = TermCounts(self.godag, self.annots)

        freq = termcounts.get_term_freq("GO:0008150")
        assert isinstance(freq, float)
        assert 0.0 <= freq <= 1.0

    def test_get_term_freq_zero_total(self):
        """Test term frequency when total count is zero."""
        # Create a mock termcounts with zero total count
        termcounts = TermCounts(self.godag, self.annots)

        # Mock the namespace to return zero total count
        with patch.object(termcounts, "get_total_count", return_value=0):
            freq = termcounts.get_term_freq("GO:0008150")
            assert freq == 0.0

    def test_get_annotations_reversed(self):
        """Test getting reversed annotations."""
        termcounts = TermCounts(self.godag, self.annots)

        reversed_annots = termcounts.get_annotations_reversed()
        assert isinstance(reversed_annots, set)
        assert len(reversed_annots) > 0

    def test_get_gosubdag_all(self):
        """Test getting GO sub-DAG with all descendants."""
        termcounts = TermCounts(self.godag, self.annots)

        gosubdag_all = termcounts.get_gosubdag_all()
        assert gosubdag_all is not None
        assert hasattr(gosubdag_all, "go2obj")

    def test_prt_objdesc(self):
        """Test printing object description."""
        termcounts = TermCounts(self.godag, self.annots)

        with patch("sys.stdout"):
            termcounts.prt_objdesc()
            # Should not raise an exception


class TestInfoContent:
    """Test cases for information content functions."""

    def setup_method(self):
        """Set up test fixtures."""
        self.godag = get_godag("go-basic.obo")
        self.annots = {
            "gene1": {"GO:0008150"},
            "gene2": {"GO:0008150"},
            "gene3": {"GO:0003674"},
        }
        self.termcounts = TermCounts(self.godag, self.annots)

    def test_get_info_content_valid(self):
        """Test getting information content for valid GO term."""
        info_content = get_info_content("GO:0008150", self.termcounts)
        assert isinstance(info_content, float)
        assert info_content >= 0.0

    def test_get_info_content_invalid(self):
        """Test getting information content for invalid GO term."""
        info_content = get_info_content("GO:9999999", self.termcounts)
        assert info_content == 0.0

    def test_get_info_content_none_termcounts(self):
        """Test getting information content with None termcounts."""
        info_content = get_info_content("GO:0008150", None)
        assert info_content == 0.0


class TestResnikSimilarity:
    """Test cases for Resnik similarity measure."""

    def setup_method(self):
        """Set up test fixtures."""
        self.godag = get_godag("go-basic.obo")
        self.annots = {
            "gene1": {"GO:0008150"},
            "gene2": {"GO:0008150"},
            "gene3": {"GO:0003674"},
        }
        self.termcounts = TermCounts(self.godag, self.annots)

    def test_resnik_sim_same_namespace(self):
        """Test Resnik similarity for terms in same namespace."""
        sim = resnik_sim("GO:0008150", "GO:0008150", self.godag, self.termcounts)
        assert isinstance(sim, float)
        assert sim >= 0.0

    def test_resnik_sim_different_namespace(self):
        """Test Resnik similarity for terms in different namespaces."""
        sim = resnik_sim("GO:0008150", "GO:0003674", self.godag, self.termcounts)
        assert sim is None

    def test_resnik_sim_invalid_terms(self):
        """Test Resnik similarity with invalid GO terms."""
        # Should handle gracefully without crashing
        try:
            resnik_sim("GO:9999999", "GO:9999998", self.godag, self.termcounts)
        except KeyError:
            # Expected for invalid GO terms
            pass


class TestLinSimilarity:
    """Test cases for Lin similarity measure."""

    def setup_method(self):
        """Set up test fixtures."""
        self.godag = get_godag("go-basic.obo")
        self.annots = {
            "gene1": {"GO:0008150"},
            "gene2": {"GO:0008150"},
            "gene3": {"GO:0003674"},
        }
        self.termcounts = TermCounts(self.godag, self.annots)

    def test_lin_sim_same_namespace(self):
        """Test Lin similarity for terms in same namespace."""
        sim = lin_sim("GO:0008150", "GO:0008150", self.godag, self.termcounts)
        assert isinstance(sim, float)
        assert 0.0 <= sim <= 1.0

    def test_lin_sim_different_namespace(self):
        """Test Lin similarity for terms in different namespaces."""
        sim = lin_sim("GO:0008150", "GO:0003674", self.godag, self.termcounts)
        assert sim is None

    def test_lin_sim_calc_with_resnik(self):
        """Test Lin similarity calculation with pre-calculated Resnik similarity."""
        sim_r = resnik_sim("GO:0008150", "GO:0008150", self.godag, self.termcounts)
        sim_l = lin_sim_calc("GO:0008150", "GO:0008150", sim_r, self.termcounts)
        assert isinstance(sim_l, float)

    def test_lin_sim_calc_same_terms(self):
        """Test Lin similarity for identical terms."""
        sim_r = resnik_sim("GO:0008150", "GO:0008150", self.godag, self.termcounts)
        sim_l = lin_sim_calc("GO:0008150", "GO:0008150", sim_r, self.termcounts)
        assert sim_l == 1.0

    def test_lin_sim_calc_zero_resnik(self):
        """Test Lin similarity with zero Resnik similarity for identical terms.

        For identical terms, Lin similarity should be 1.0 even when Resnik similarity is 0.0,
        because identical terms are maximally similar by definition.
        """
        sim_l = lin_sim_calc("GO:0008150", "GO:0008150", 0.0, self.termcounts)
        assert sim_l == 1.0

    def test_lin_sim_calc_with_default(self):
        """Test Lin similarity with default value."""
        sim_l = lin_sim_calc(
            "GO:0008150", "GO:0003674", None, self.termcounts, dfltval=0.5
        )
        assert sim_l == 0.5


class TestSchlickerSimilarity:
    """Test cases for Schlicker similarity measure."""

    def setup_method(self):
        """Set up test fixtures."""
        self.godag = get_godag("go-basic.obo")
        self.annots = {
            "gene1": {"GO:0008150"},
            "gene2": {"GO:0008150"},
            "gene3": {"GO:0003674"},
        }
        self.termcounts = TermCounts(self.godag, self.annots)

    def test_schlicker_sim_same_namespace(self):
        """Test Schlicker similarity for terms in same namespace."""
        sim = schlicker_sim("GO:0008150", "GO:0008150", self.godag, self.termcounts)
        assert isinstance(sim, float)
        assert 0.0 <= sim <= 1.0

    def test_schlicker_sim_different_namespace(self):
        """Test Schlicker similarity for terms in different namespaces."""
        sim = schlicker_sim("GO:0008150", "GO:0003674", self.godag, self.termcounts)
        assert sim is None

    def test_schlicker_sim_calc_with_components(self):
        """Test Schlicker similarity calculation with pre-calculated components."""
        sim_r = resnik_sim("GO:0008150", "GO:0008150", self.godag, self.termcounts)
        tfreq = get_freq_msca("GO:0008150", "GO:0008150", self.godag, self.termcounts)
        sim_s = schlicker_sim_calc(
            "GO:0008150", "GO:0008150", sim_r, tfreq, self.termcounts
        )
        assert isinstance(sim_s, float)

    def test_schlicker_sim_calc_same_terms(self):
        """Test Schlicker similarity for identical terms."""
        sim_r = resnik_sim("GO:0008150", "GO:0008150", self.godag, self.termcounts)
        tfreq = get_freq_msca("GO:0008150", "GO:0008150", self.godag, self.termcounts)
        sim_s = schlicker_sim_calc(
            "GO:0008150", "GO:0008150", sim_r, tfreq, self.termcounts
        )
        assert sim_s == 1.0 - tfreq

    def test_schlicker_sim_calc_zero_resnik(self):
        """Test Schlicker similarity with zero Resnik similarity."""
        tfreq = get_freq_msca("GO:0008150", "GO:0008150", self.godag, self.termcounts)
        sim_s = schlicker_sim_calc(
            "GO:0008150", "GO:0008150", 0.0, tfreq, self.termcounts
        )
        assert sim_s == 0.0


class TestFrequencyFunctions:
    """Test cases for frequency-related functions."""

    def setup_method(self):
        """Set up test fixtures."""
        self.godag = get_godag("go-basic.obo")
        self.annots = {
            "gene1": {"GO:0008150"},
            "gene2": {"GO:0008150"},
            "gene3": {"GO:0003674"},
        }
        self.termcounts = TermCounts(self.godag, self.annots)

    def test_get_freq_msca_same_namespace(self):
        """Test getting frequency of MSCA for terms in same namespace."""
        freq = get_freq_msca("GO:0008150", "GO:0008150", self.godag, self.termcounts)
        assert isinstance(freq, float)
        assert 0.0 <= freq <= 1.0

    def test_get_freq_msca_different_namespace(self):
        """Test getting frequency of MSCA for terms in different namespaces."""
        freq = get_freq_msca("GO:0008150", "GO:0003674", self.godag, self.termcounts)
        assert freq == 0

    def test_get_freq_msca_invalid_terms(self):
        """Test getting frequency of MSCA for invalid terms."""
        freq = get_freq_msca("GO:9999999", "GO:9999998", self.godag, self.termcounts)
        assert freq == 0


class TestAncestorFunctions:
    """Test cases for ancestor-related functions."""

    def setup_method(self):
        """Set up test fixtures."""
        self.godag = get_godag("go-basic.obo")

    def test_common_parent_go_ids_single(self):
        """Test common parent GO IDs with single term."""
        goids = ["GO:0008150"]
        common_parents = common_parent_go_ids(goids, self.godag)
        assert isinstance(common_parents, set)
        assert len(common_parents) > 0

    def test_common_parent_go_ids_multiple(self):
        """Test common parent GO IDs with multiple terms."""
        goids = ["GO:0008150", "GO:0008150"]
        common_parents = common_parent_go_ids(goids, self.godag)
        assert isinstance(common_parents, set)
        assert len(common_parents) > 0

    def test_common_parent_go_ids_invalid(self):
        """Test common parent GO IDs with invalid terms."""
        goids = ["GO:9999999", "GO:9999998"]
        # Should handle gracefully without crashing
        try:
            common_parents = common_parent_go_ids(goids, self.godag)
            assert isinstance(common_parents, set)
        except KeyError:
            # This is expected for invalid GO terms
            pass

    def test_deepest_common_ancestor_single(self):
        """Test deepest common ancestor with single term."""
        goids = ["GO:0008150"]
        dca = deepest_common_ancestor(goids, self.godag)
        assert dca == "GO:0008150"

    def test_deepest_common_ancestor_multiple(self):
        """Test deepest common ancestor with multiple terms."""
        goids = ["GO:0008150", "GO:0008150"]
        dca = deepest_common_ancestor(goids, self.godag)
        assert isinstance(dca, str)
        assert dca.startswith("GO:")

    def test_deepest_common_ancestor_invalid(self):
        """Test deepest common ancestor with invalid terms."""
        goids = ["GO:9999999", "GO:9999998"]
        # Should handle gracefully without crashing
        try:
            dca = deepest_common_ancestor(goids, self.godag)
            assert isinstance(dca, str)
        except (KeyError, ValueError):
            # This is expected for invalid GO terms
            pass


class TestDistanceFunctions:
    """Test cases for distance and similarity functions."""

    def setup_method(self):
        """Set up test fixtures."""
        self.godag = get_godag("go-basic.obo")

    def test_min_branch_length_same_namespace(self):
        """Test minimum branch length for terms in same namespace."""
        length = min_branch_length("GO:0008150", "GO:0008150", self.godag, None)
        assert isinstance(length, int)
        assert length >= 0

    def test_min_branch_length_different_namespace(self):
        """Test minimum branch length for terms in different namespaces."""
        length = min_branch_length("GO:0008150", "GO:0003674", self.godag, None)
        assert length is None

    def test_min_branch_length_with_branch_dist(self):
        """Test minimum branch length with branch distance."""
        branch_dist = 2
        length = min_branch_length("GO:0008150", "GO:0003674", self.godag, branch_dist)
        assert isinstance(length, int)
        assert length > 0

    def test_semantic_distance_same_namespace(self):
        """Test semantic distance for terms in same namespace."""
        distance = semantic_distance("GO:0008150", "GO:0008150", self.godag)
        assert isinstance(distance, int)
        assert distance >= 0

    def test_semantic_distance_different_namespace(self):
        """Test semantic distance for terms in different namespaces."""
        distance = semantic_distance("GO:0008150", "GO:0003674", self.godag)
        assert distance is None

    def test_semantic_distance_with_branch_dist(self):
        """Test semantic distance with branch distance."""
        branch_dist = 3
        distance = semantic_distance(
            "GO:0008150", "GO:0003674", self.godag, branch_dist
        )
        assert isinstance(distance, int)
        assert distance > 0

    def test_semantic_similarity_same_namespace(self):
        """Test semantic similarity for terms in same namespace."""
        similarity = semantic_similarity("GO:0008150", "GO:0008150", self.godag)
        assert isinstance(similarity, float)
        assert similarity == 1.0  # Same term should have similarity 1.0

    def test_semantic_similarity_different_namespace(self):
        """Test semantic similarity for terms in different namespaces."""
        similarity = semantic_similarity("GO:0008150", "GO:0003674", self.godag)
        assert similarity is None

    def test_semantic_similarity_with_branch_dist(self):
        """Test semantic similarity with branch distance."""
        branch_dist = 2
        similarity = semantic_similarity(
            "GO:0008150", "GO:0003674", self.godag, branch_dist
        )
        assert isinstance(similarity, float)
        assert 0.0 < similarity <= 1.0

    def test_semantic_similarity_zero_distance(self):
        """Test semantic similarity when distance is zero."""
        # Same term should have distance 0 and similarity 1.0
        similarity = semantic_similarity("GO:0008150", "GO:0008150", self.godag)
        assert similarity == 1.0


class TestEdgeCases:
    """Test cases for edge cases and error conditions."""

    def setup_method(self):
        """Set up test fixtures."""
        self.godag = get_godag("go-basic.obo")
        self.annots = {
            "gene1": {"GO:0008150"},
            "gene2": {"GO:0003674"},
        }
        self.termcounts = TermCounts(self.godag, self.annots)

    def test_empty_annotations(self):
        """Test with empty annotations."""
        empty_annots = {}
        termcounts = TermCounts(self.godag, empty_annots)
        assert len(termcounts.gocnts) == 0
        assert len(termcounts.gene2gos) == 0

    def test_none_annotations(self):
        """Test with None annotations."""
        with pytest.raises((TypeError, AttributeError)):
            TermCounts(self.godag, None)

    def test_invalid_godag(self):
        """Test with invalid GODag."""
        with pytest.raises((TypeError, AttributeError)):
            TermCounts(None, self.annots)

    def test_mixed_valid_invalid_annotations(self):
        """Test with mix of valid and invalid GO terms in annotations."""
        mixed_annots = {
            "gene1": {"GO:0008150", "GO:9999999"},  # Valid and invalid GO term
            "gene2": {"GO:0003674"},
        }
        # Should handle gracefully, ignoring invalid terms
        termcounts = TermCounts(self.godag, mixed_annots)
        assert len(termcounts.gocnts) > 0

    def test_very_large_annotations(self):
        """Test with very large annotation set."""
        large_annots = {}
        for i in range(1000):
            large_annots[f"gene{i}"] = {"GO:0008150", "GO:0003674"}

        # Should handle large datasets without issues
        termcounts = TermCounts(self.godag, large_annots)
        assert len(termcounts.gocnts) > 0
        assert len(termcounts.gene2gos) == 1000


class TestIntegration:
    """Integration tests using real GO data."""

    def setup_method(self):
        """Set up test fixtures."""
        self.godag = get_godag("go-basic.obo")

    def test_real_go_terms_bp(self):
        """Test with real biological process GO terms."""
        # Use real GO terms from the ontology
        bp_terms = [
            "GO:0008150",  # biological_process
            "GO:0009987",  # cellular process
            "GO:0016043",  # cellular component organization
        ]

        # Create annotations with these terms
        annots = {
            "gene1": {bp_terms[0]},
            "gene2": {bp_terms[1]},
            "gene3": {bp_terms[2]},
        }

        termcounts = TermCounts(self.godag, annots)

        # Test similarity between related terms
        sim = semantic_similarity(bp_terms[1], bp_terms[2], self.godag)
        assert isinstance(sim, float)
        assert 0.0 <= sim <= 1.0

    def test_real_go_terms_cross_namespace(self):
        """Test similarity across different namespaces."""
        bp_term = "GO:0008150"  # biological_process
        mf_term = "GO:0003674"  # molecular_function
        cc_term = "GO:0005575"  # cellular_component

        annots = {
            "gene1": {bp_term},
            "gene2": {mf_term},
            "gene3": {cc_term},
        }

        termcounts = TermCounts(self.godag, annots)

        # Cross-namespace similarities should be None
        sim_bp_mf = semantic_similarity(bp_term, mf_term, self.godag)
        sim_bp_cc = semantic_similarity(bp_term, cc_term, self.godag)
        sim_mf_cc = semantic_similarity(mf_term, cc_term, self.godag)

        assert sim_bp_mf is None
        assert sim_bp_cc is None
        assert sim_mf_cc is None

    def test_term_frequency_calculation(self):
        """Test term frequency calculation with real data."""
        annots = {
            "gene1": {"GO:0008150"},
            "gene2": {"GO:0008150"},
            "gene3": {"GO:0003674"},
            "gene4": {"GO:0003674"},
            "gene5": {"GO:0005575"},
        }

        termcounts = TermCounts(self.godag, annots)

        # Test frequency calculation
        freq_bp = termcounts.get_term_freq("GO:0008150")
        freq_mf = termcounts.get_term_freq("GO:0003674")
        freq_cc = termcounts.get_term_freq("GO:0005575")

        assert isinstance(freq_bp, float)
        assert isinstance(freq_mf, float)
        assert isinstance(freq_cc, float)
        assert 0.0 <= freq_bp <= 1.0
        assert 0.0 <= freq_mf <= 1.0
        assert 0.0 <= freq_cc <= 1.0

    def test_comprehensive_similarity_measures(self):
        """Test all similarity measures with real data."""
        annots = {
            "gene1": {"GO:0008150"},
            "gene2": {"GO:0008150"},
            "gene3": {"GO:0009987"},
        }

        termcounts = TermCounts(self.godag, annots)

        go1, go2 = "GO:0008150", "GO:0009987"

        # Test all similarity measures
        resnik = resnik_sim(go1, go2, self.godag, termcounts)
        lin = lin_sim(go1, go2, self.godag, termcounts)
        schlicker = schlicker_sim(go1, go2, self.godag, termcounts)

        if resnik is not None:
            assert isinstance(resnik, float)
            assert resnik >= 0.0

        if lin is not None:
            assert isinstance(lin, float)
            assert 0.0 <= lin <= 1.0

        if schlicker is not None:
            assert isinstance(schlicker, float)
            assert 0.0 <= schlicker <= 1.0


if __name__ == "__main__":
    # Run the tests
    pytest.main([__file__, "-v"])
