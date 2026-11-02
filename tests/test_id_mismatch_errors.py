"""IDs that don't match the annotations fail with a message that shows the mismatch (#254, #260, #186)."""

import os
import subprocess
import sys

import pytest

from goatools.anno.genetogo_reader import Gene2GoReader

MINI_OBO = os.path.join(os.path.dirname(os.path.abspath(__file__)), "data", "mini_obo.obo")


def _run(tmp_path, study):
    """Run find_enrichment with gene symbols against UniProt-style annotation IDs."""
    (tmp_path / "assoc.txt").write_text("P11\tGO:0000004\nP12\tGO:0000005\nP13\tGO:0000006\n")
    (tmp_path / "pop.txt").write_text("geneA\ngeneB\ngeneC\ngeneD\n")
    (tmp_path / "study.txt").write_text("\n".join(study) + "\n")
    cmd = [sys.executable, "-m", "goatools", "find_enrichment"] + [
        str(tmp_path / f) for f in ("study.txt", "pop.txt", "assoc.txt")] + [
        "--obo=" + MINI_OBO, "--outfile=" + str(tmp_path / "out.tsv")]
    res = subprocess.run(cmd, capture_output=True, text=True, env=dict(os.environ, COLUMNS="300"))
    return res.returncode, " ".join((res.stdout + res.stderr).split())


def test_population_not_in_annotations(tmp_path):
    """Used to log a fatal error and exit 0."""
    code, out = _run(tmp_path, ["geneA", "geneB"])
    assert code == 1
    assert "Population IDs (e.g. geneA, geneB, geneC) must match the annotation IDs (e.g. P11, P12, P13)" in out


def test_study_not_in_population(tmp_path):
    """Names --min_overlap and the missing study IDs."""
    code, out = _run(tmp_path, ["x1", "x2", "x3", "geneA"])
    assert code == 1
    assert "only 25% of genes/proteins in the study" in out
    assert "(--min_overlap 0.7). Study IDs not in the population include: x1, x2, x3" in out


def test_gene2go_gene_ids_as_taxids(tmp_path):
    """#186: gene IDs passed as taxids= named, with what taxids are."""
    fin = tmp_path / "gene2go"
    fin.write_text("#tax_id\tGeneID\tGO_ID\tEvidence\tQualifier\tGO_term\tPubMed\tCategory\n"
                   "9606\t1634\tGO:0005515\tIPI\tenables\tprotein binding\t12345\tFunction\n")
    assert list(Gene2GoReader(str(fin), taxids=[9606]).taxid2asscs) == [9606]
    with pytest.raises(AssertionError, match=r"taxids\(\[1634, 347\]\).*not gene IDs"):
        Gene2GoReader(str(fin), taxids=[1634, 347])
