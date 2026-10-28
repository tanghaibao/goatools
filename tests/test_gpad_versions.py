"""GPAD 2.0 files read to the same annotations as the equivalent GPAD 1.1 files.

GPAD 2.0 (the format of the GO Consortium's current goa_human.gpad) keeps the 12
columns but joins DB:DB_Object_ID, adds a negation column, replaces the Qualifier
label with a relation ID, and writes dates as YYYY-MM-DD.
"""

__copyright__ = "Copyright (C) 2016-present, DV Klopfenstein, H Tang. All rights reserved."

from datetime import date

import pytest

from goatools.anno.gpad_reader import GpadReader

# The same four annotations in each format
GPAD_1_1 = """!gpa-version: 1.1
UniProtKB\tP1\tenables\tGO:0005515\tPMID:1\tECO:0000353\tUniProtKB:Q1\t\t20260425\tIntAct\t\t
UniProtKB\tP2\tNOT|involved_in\tGO:0006915\tPMID:2|GO_REF:0000024\tECO:0000250\t\t\t20260101\tUniProt\t\t
MGI\tMGI:101757\tlocated_in\tGO:0005737\tPMID:3\tECO:0000314\t\ttaxon:10090\t20260102\tMGI\t\t
ComplexPortal\tCPX-1\tpart_of\tGO:0005634\tGO_REF:0000107\tECO:0000265\t\t\t20260103\tComplexPortal\t\t
"""

GPAD_2_0 = """!gpad-version: 2.0
!generated-by: GOC
UniProtKB:P1\t\tRO:0002327\tGO:0005515\tPMID:1\tECO:0000353\tUniProtKB:Q1\t\t2026-04-25\tIntAct\t\t
UniProtKB:P2\tNOT\tRO:0002331\tGO:0006915\tPMID:2|GO_REF:0000024\tECO:0000250\t\t\t2026-01-01\tUniProt\t\t
MGI:MGI:101757\t\tRO:0001025\tGO:0005737\tPMID:3\tECO:0000314\t\tNCBITaxon:10090\t2026-01-02\tMGI\t\t
ComplexPortal:CPX-1\t\tBFO:0000050\tGO:0005634\tGO_REF:0000107\tECO:0000265\t\t\t2026-01-03\tComplexPortal\t\t
"""


def _read(tmp_path, name, text):
    fin = tmp_path / name
    fin.write_text(text)
    return GpadReader(str(fin))


def test_gpad_2_0_matches_1_1(tmp_path):
    """Every field of every annotation agrees across the two formats."""
    v11 = _read(tmp_path, "v11.gpa", GPAD_1_1).associations
    v20 = _read(tmp_path, "v20.gpad", GPAD_2_0).associations
    # every annotation is read, including the first one after the header
    assert [nt.DB_ID for nt in v11] == ["P1", "P2", "MGI:101757", "CPX-1"]
    assert v20 == v11


def test_gpad_2_0_fields(tmp_path):
    """The 2.0 columns that changed are converted, not just passed through."""
    p1, p2, mgi, cpx = _read(tmp_path, "v20.gpad", GPAD_2_0).associations
    assert (p1.DB, p1.DB_ID, p1.Qualifier) == ("UniProtKB", "P1", {"enables"})
    assert p2.Qualifier == {"NOT", "involved_in"}  # negation column
    assert p2.DB_Reference == {"PMID:2", "GO_REF:0000024"}
    assert (mgi.DB, mgi.DB_ID, mgi.Taxon) == ("MGI", "MGI:101757", 10090)
    assert mgi.Date == date(2026, 1, 2)
    assert cpx.Qualifier == {"part_of"}


def test_gpad_2_0_taxa_and_properties(tmp_path):
    """2.0 allows several interacting taxa and repeated annotation properties."""
    props = ("contributor-id=orcid:0000-0002-1478-7671|model-state=production|"
             "contributor-id=orcid:0000-0001-7476-6306")
    text = GPAD_2_0.replace(
        "\t\t2026-04-25\tIntAct\t\t\n",
        "\tNCBITaxon:9606|NCBITaxon:10090\t2026-04-25\tIntAct\t\t" + props + "\n")
    p1, p2, mgi, _ = _read(tmp_path, "v20.gpad", text).associations
    assert p1.Taxon == [9606, 10090]
    assert p1.Properties == {
        "contributor-id": ["orcid:0000-0002-1478-7671", "orcid:0000-0001-7476-6306"],
        "model-state": ["production"],
    }
    # one taxon stays an int and no properties stay empty, as in 1.1
    assert (mgi.Taxon, mgi.Properties, p2.Taxon) == (10090, {}, None)


def test_gpad_2_0_id2gos(tmp_path):
    """Downstream id2gos is unchanged; the NOT annotation is excluded as in 1.1."""
    id2gos_11 = _read(tmp_path, "v11.gpa", GPAD_1_1).get_id2gos(prt=None)
    id2gos_20 = _read(tmp_path, "v20.gpad", GPAD_2_0).get_id2gos(prt=None)
    assert id2gos_20 == id2gos_11
    assert "P2" not in id2gos_20 and "P1" in id2gos_20


def test_gpad_2_0_unknown_relation(tmp_path, capsys):
    """A relation with no GPAD 1.1 qualifier is reported, not silently dropped."""
    text = GPAD_2_0.replace("RO:0002327", "RO:9999999")
    with pytest.raises(SystemExit):
        _read(tmp_path, "bad.gpad", text)
    assert "UNRECOGNIZED RELATION(RO:9999999)" in capsys.readouterr().out
