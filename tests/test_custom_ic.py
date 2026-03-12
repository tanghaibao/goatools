#!/usr/bin/env python
"""Test information content calculation for custom GOA in id2go format.

Verifies that users can compute per-term information content (IC) from
a custom id2go annotation file (or an in-memory dict) for any species,
as documented in notebooks/custom_tinfo.ipynb.
"""

import os
import tempfile

import pytest

from goatools.anno.idtogos_reader import IdToGosReader
from goatools.base import get_godag as base_get_godag
from goatools.semantic import TermCounts, get_info_content, lin_sim, resnik_sim

# Use the small heartjogging OBO file that ships with the test suite so that
# this test never requires a network connection.
_DIR_TESTS = os.path.dirname(os.path.abspath(__file__))
_REPO = os.path.abspath(os.path.join(_DIR_TESTS, ".."))
_OBO = os.path.join(_REPO, "tests", "data", "heartjogging.obo")


def _get_godag():
    """Load the heartjogging test ontology (no network required)."""
    godag_raw = base_get_godag(_OBO)
    return {o.item_id: o for o in godag_raw.values()}


# ---------------------------------------------------------------------------
# Custom annotations used throughout the tests
# Two 'tools' annotated the same gene-set differently; we compare their ICs.
# ---------------------------------------------------------------------------
# Annotations from a first annotation tool (more specific terms preferred)
_ANNOTS_TOOL1 = {
    "geneA": {"GO:0007507", "GO:0007368"},  # heart development, det. L/R symmetry
    "geneB": {"GO:0007507", "GO:0007275"},  # heart development, multicell. dev.
    "geneC": {"GO:0003007", "GO:0007275"},  # heart morphogenesis, multicell. dev.
}

# Annotations from a second annotation tool (uses broader/different terms)
_ANNOTS_TOOL2 = {
    "geneA": {"GO:0007275"},               # multicellular organism development only
    "geneB": {"GO:0007507"},               # heart development only
    "geneC": {"GO:0007275"},               # multicellular organism development only
}


class TestCustomICFromDict:
    """IC calculation from an in-memory dict (no file I/O needed)."""

    def setup_method(self):
        self.godag = _get_godag()
        self.tcnt1 = TermCounts(self.godag, _ANNOTS_TOOL1)
        self.tcnt2 = TermCounts(self.godag, _ANNOTS_TOOL2)

    def test_termcounts_created(self):
        """TermCounts can be built from a plain gene->GO-set dict."""
        assert len(self.tcnt1.gocnts) > 0
        assert len(self.tcnt1.gene2gos) == 3

    def test_ic_non_negative(self):
        """IC values must be >= 0 for all annotated terms."""
        for go_id in _ANNOTS_TOOL1:
            # go_id here is a gene; iterate over GO IDs instead
            pass
        for gene, goids in _ANNOTS_TOOL1.items():
            for go_id in goids:
                ic = get_info_content(go_id, self.tcnt1)
                assert ic >= 0.0, f"Negative IC for {go_id}"

    def test_ic_more_specific_terms_higher(self):
        """Less frequently annotated terms have higher IC than commonly annotated ones.

        In _ANNOTS_TOOL1:
          GO:0007368 (determination of L/R symmetry) is annotated to 1 gene (geneA)
          GO:0007275 (multicellular organism development) is annotated to 2 genes (geneB, geneC)

        frequency(GO:0007368) = 1/3 < frequency(GO:0007275) = 2/3
        => IC(GO:0007368) > IC(GO:0007275)
        """
        ic_rare = get_info_content("GO:0007368", self.tcnt1)    # 1 gene → rare
        ic_common = get_info_content("GO:0007275", self.tcnt1)  # 2 genes → common
        assert ic_rare > ic_common, (
            f"Expected IC(det. L/R symmetry)={ic_rare:.4f} > "
            f"IC(multicell dev)={ic_common:.4f}"
        )

    def test_ic_unannotated_term_is_zero(self):
        """A GO term absent from any annotation should return IC=0."""
        ic = get_info_content("GO:0008150", self.tcnt1)  # root BP term
        # Root term gets propagated count == total, so tfreq == 1 => -log(1) == 0
        assert ic == 0.0

    def test_ic_differs_between_annotation_sets(self):
        """IC values depend on the annotation set; two tools give different ICs."""
        ic1 = get_info_content("GO:0007507", self.tcnt1)
        ic2 = get_info_content("GO:0007507", self.tcnt2)
        # tool2 has only one gene annotated to heart development out of 3 total
        # tool1 has two genes -> different frequencies
        assert ic1 != ic2, "Expected different IC values for different annotation sets"

    def test_resnik_sim_custom(self):
        """Resnik similarity can be calculated with custom annotations."""
        go1, go2 = "GO:0007507", "GO:0003007"
        sim = resnik_sim(go1, go2, self.godag, self.tcnt1)
        assert sim is not None
        assert sim >= 0.0

    def test_lin_sim_custom(self):
        """Lin similarity can be calculated with custom annotations."""
        go1, go2 = "GO:0007507", "GO:0003007"
        sim = lin_sim(go1, go2, self.godag, self.tcnt1)
        # Both terms are annotated so Lin sim should be defined and in [0, 1]
        assert sim is not None
        assert 0.0 <= sim <= 1.0

    def test_compare_two_annotation_sets(self):
        """Helper: compare IC for shared GO terms across two annotation sets."""
        shared_goids = set.union(*_ANNOTS_TOOL1.values()) & set.union(
            *_ANNOTS_TOOL2.values()
        )
        assert shared_goids, "Expected at least one shared GO term"
        for go_id in shared_goids:
            ic1 = get_info_content(go_id, self.tcnt1)
            ic2 = get_info_content(go_id, self.tcnt2)
            # ICs must be valid floats
            assert isinstance(ic1, float)
            assert isinstance(ic2, float)


class TestCustomICFromIdToGosFile:
    """IC calculation by reading a custom id2go annotation file."""

    def setup_method(self):
        self.godag = _get_godag()

    def _write_and_read(self, id2gos):
        """Write id2gos to a temp file and read back via IdToGosReader."""
        with tempfile.NamedTemporaryFile(
            mode="w", suffix=".txt", delete=False
        ) as fh:
            tmp_path = fh.name
        try:
            IdToGosReader.wr_id2gos(tmp_path, id2gos)
            reader = IdToGosReader(tmp_path, godag=self.godag)
            return reader.get_id2gos()
        finally:
            os.unlink(tmp_path)

    def test_roundtrip_preserves_annotations(self):
        """Writing then reading an id2go file returns the original gene->GO mapping."""
        loaded = self._write_and_read(_ANNOTS_TOOL1)
        for gene, goids in _ANNOTS_TOOL1.items():
            assert gene in loaded, f"Gene {gene} missing after roundtrip"
            assert goids == loaded[gene], (
                f"GO set mismatch for {gene}: expected {goids}, got {loaded[gene]}"
            )

    def test_termcounts_from_file(self):
        """TermCounts built from a file-loaded id2go annotation is valid."""
        loaded = self._write_and_read(_ANNOTS_TOOL1)
        tcnt = TermCounts(self.godag, loaded)
        assert len(tcnt.gocnts) > 0
        # heart development should have count == 2 (geneA and geneB annotated it)
        assert tcnt.get_count("GO:0007507") == 2

    def test_ic_from_file_matches_dict(self):
        """IC from a file-loaded annotation must equal IC from the in-memory dict."""
        loaded = self._write_and_read(_ANNOTS_TOOL1)
        tcnt_file = TermCounts(self.godag, loaded)
        tcnt_dict = TermCounts(self.godag, _ANNOTS_TOOL1)

        for gene, goids in _ANNOTS_TOOL1.items():
            for go_id in goids:
                ic_file = get_info_content(go_id, tcnt_file)
                ic_dict = get_info_content(go_id, tcnt_dict)
                assert ic_file == pytest.approx(ic_dict), (
                    f"IC mismatch for {go_id}: file={ic_file}, dict={ic_dict}"
                )

    def test_compare_two_custom_goas_from_files(self):
        """IC comparison between two custom GOAs loaded from files."""
        loaded1 = self._write_and_read(_ANNOTS_TOOL1)
        loaded2 = self._write_and_read(_ANNOTS_TOOL2)

        tcnt1 = TermCounts(self.godag, loaded1)
        tcnt2 = TermCounts(self.godag, loaded2)

        # GO:0007507 (heart development): tool1 has 2/3 genes, tool2 has 1/3
        ic1 = get_info_content("GO:0007507", tcnt1)
        ic2 = get_info_content("GO:0007507", tcnt2)

        # tool1 annotates 2 genes -> higher frequency -> lower IC
        assert ic1 < ic2, (
            f"Expected tool1 IC ({ic1:.4f}) < tool2 IC ({ic2:.4f}) for GO:0007507"
        )
