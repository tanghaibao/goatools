"""Regression tests for indexing dict.items() views (a TypeError on Python 3)."""

import pytest

from goatools.parsers.ncbi_gene_file_reader import NCBIgeneFileReader
from goatools.ratio import _chk_gene2go


def test_ncbi_gene_file_reader_chromosome_column(tmp_path):
    """Values in a 'chromosome' column are right-aligned for sorting."""
    fin = tmp_path / "genes.tsv"
    fin.write_text("GeneID\tSymbol\tchromosome\n1017\tCDK2\t9\n7157\tTP53\t17\n")
    data, _ = NCBIgeneFileReader(str(fin), log=None).run("lists", None)
    assert data == [[1017, "CDK2", " 9"], [7157, "TP53", "17"]]


def test_chk_gene2go_rejects_go2gene():
    """go2gene associations raise the intended error, not a TypeError."""
    with pytest.raises(Exception, match="EXPECTED TO BE gene2go") as excinfo:
        _chk_gene2go({"GO:0008150": {"CDK2"}})
    assert not isinstance(excinfo.value, TypeError)
