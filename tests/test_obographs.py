"""GODag reads OBO Graphs JSON (go-basic.json) into the same DAG as the equivalent obo file."""

import json

import pytest

from goatools.obo_parser import GODag

OPTIONAL_ATTRS = {"def", "synonym", "relationship", "xref", "subset", "comment", "consider", "replaced_by"}

OBO = """format-version: 1.2
data-version: releases/2026-07-26
default-namespace: gene_ontology

[Term]
id: GO:0000001
name: mitochondrion inheritance
namespace: biological_process
alt_id: GO:0000099
def: "The distribution of mitochondria into daughter cells." [GOC:mcc, PMID:10873824]
comment: A comment.
subset: goslim_generic
synonym: "Mito inheritance" RELATED systematic_synonym [GOC:x]
synonym: "mitochondrial inheritance" EXACT []
xref: Wikipedia:Mitochondrion
is_a: GO:0008150 ! biological_process
relationship: part_of GO:0000002 ! mitochondrial genome maintenance

[Term]
id: GO:0000002
name: mitochondrial genome maintenance
namespace: biological_process
is_a: GO:0008150 ! biological_process

[Term]
id: GO:0000003
name: obsolete reproduction
namespace: biological_process
is_obsolete: true
replaced_by: GO:0008150
consider: GO:0000002

[Term]
id: GO:0008150
name: biological_process
namespace: biological_process

[Typedef]
id: part_of
name: part of
namespace: external
xref: BFO:0000050
"""

P = "http://purl.obolibrary.org/obo/"
O = "http://www.geneontology.org/formats/oboInOwl#"
BP = {"pred": O + "hasOBONamespace", "val": "biological_process"}

# The same content as OBO, in the order and form go-basic.json uses
GRAPH = {"graphs": [{
    "id": P + "go.owl",
    "meta": {
        "version": P + "go/releases/2026-07-26/go.owl",
        "basicPropertyValues": [
            {"pred": O + "hasOBOFormatVersion", "val": "1.2"},
            {"pred": O + "default-namespace", "val": "gene_ontology"},
        ],
    },
    "nodes": [
        {"id": P + "GO_0008150", "lbl": "biological_process", "type": "CLASS",
         "meta": {"basicPropertyValues": [BP]}},
        {"id": P + "GO_0000001", "lbl": "mitochondrion inheritance", "type": "CLASS", "meta": {
            "definition": {"val": "The distribution of mitochondria into daughter cells.",
                           "xrefs": ["PMID:10873824", "GOC:mcc"]},
            "comments": ["A comment."],
            "subsets": [P + "go#goslim_generic"],
            "synonyms": [
                {"pred": "hasExactSynonym", "val": "mitochondrial inheritance"},
                {"pred": "hasRelatedSynonym", "val": "Mito inheritance", "xrefs": ["GOC:x"],
                 "synonymType": P + "go#systematic_synonym"},
            ],
            "xrefs": [{"val": "Wikipedia:Mitochondrion"}],
            "basicPropertyValues": [BP, {"pred": O + "hasAlternativeId", "val": "GO:0000099"}],
        }},
        # An alt_id is also a stub class; it is not a term of its own
        {"id": P + "GO_0000099", "type": "CLASS", "meta": {
            "deprecated": True,
            "basicPropertyValues": [{"pred": P + "IAO_0100001", "val": P + "GO_0000001"}],
        }},
        {"id": P + "GO_0000002", "lbl": "mitochondrial genome maintenance", "type": "CLASS",
         "meta": {"basicPropertyValues": [BP]}},
        {"id": P + "GO_0000003", "lbl": "obsolete reproduction", "type": "CLASS", "meta": {
            "deprecated": True,
            "basicPropertyValues": [
                {"pred": P + "IAO_0100001", "val": P + "GO_0008150"},
                BP,
                {"pred": O + "consider", "val": "GO:0000002"},
            ],
        }},
        {"id": P + "BFO_0000050", "lbl": "part of", "type": "PROPERTY", "meta": {
            "basicPropertyValues": [
                {"pred": O + "hasOBONamespace", "val": "external"},
                {"pred": O + "shorthand", "val": "part_of"},
            ],
        }},
    ],
    "edges": [
        {"sub": P + "GO_0000001", "pred": P + "BFO_0000050", "obj": P + "GO_0000002"},
        {"sub": P + "GO_0000001", "pred": "is_a", "obj": P + "GO_0008150"},
        {"sub": P + "GO_0000002", "pred": "is_a", "obj": P + "GO_0008150"},
    ],
}]}


def _summary(rec):
    """Every GOTerm field, with linked records as GO IDs."""
    vals = dict(vars(rec))
    for attr in ("parents", "children"):
        vals[attr] = sorted(r.item_id for r in vals[attr])
    for attr in ("relationship", "relationship_rev"):
        vals[attr] = {k: sorted(r.item_id for r in v) for k, v in vals[attr].items()}
    return vals


@pytest.fixture(name="dags")
def _dags(tmp_path):
    """(obo GODag, json GODag) for each load_obsolete setting."""
    fobo, fjson = tmp_path / "go.obo", tmp_path / "go.json"
    fobo.write_text(OBO)
    fjson.write_text(json.dumps(GRAPH))
    return {
        obsolete: tuple(
            GODag(str(fin), optional_attrs=OPTIONAL_ATTRS, load_obsolete=obsolete, prt=None)
            for fin in (fobo, fjson))
        for obsolete in (False, True)
    }


@pytest.mark.parametrize("load_obsolete", [False, True])
def test_json_matches_obo(dags, load_obsolete):
    """Same terms, alt IDs, fields, typedefs and release, from either format."""
    dag_obo, dag_json = dags[load_obsolete]
    assert list(dag_json) == list(dag_obo)
    for goid, rec in dag_obo.items():
        assert _summary(dag_json[goid]) == _summary(rec), goid
    assert {k: vars(v) for k, v in dag_json.typedefs.items()} == \
        {k: vars(v) for k, v in dag_obo.typedefs.items()}
    assert dag_json.data_version == dag_obo.data_version == "releases/2026-07-26"
    assert dag_json.version.split(": ", 1)[1] == dag_obo.version.split(": ", 1)[1]


def test_json_fields(dags):
    """Spot-check the conversions that differ between the formats."""
    dag = dags[True][1]
    term = dag["GO:0000001"]
    assert dag["GO:0000099"] is term
    assert [s.text for s in term.synonym] == ["Mito inheritance", "mitochondrial inheritance"]
    assert term.synonym[0][1:3] == ("RELATED", "systematic_synonym")
    assert term.defn == '"The distribution of mitochondria into daughter cells." [GOC:mcc, PMID:10873824]'
    assert term.relationship == {"part_of": {dag["GO:0000002"]}}
    assert (term.subset, term.xref) == ({"goslim_generic"}, {"Wikipedia:Mitochondrion"})
    obsolete = dag["GO:0000003"]
    assert obsolete.is_obsolete
    assert (obsolete.replaced_by, obsolete.consider) == ("GO:0008150", {"GO:0000002"})
