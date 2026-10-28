"""Read a GO DAG from OBO Graphs JSON (go-basic.json) into the same records as an obo file."""

__copyright__ = "Copyright (C) 2010-present, H Tang et al. All rights reserved."

import collections as cx
import json
import re
from typing import Dict, Iterator, List, Set

from goatools.godag.typedef import TypeDef, add_to_typedef
from goatools.obo_parser import GOTerm, OBOReader

OBO_PURL: str = "http://purl.obolibrary.org/obo/"
OBOINOWL: str = "http://www.geneontology.org/formats/oboInOwl#"

# Node property predicate -> obo tag
PRED2TAG: Dict[str, str] = {
    OBOINOWL + "hasOBONamespace": "namespace",
    OBOINOWL + "hasAlternativeId": "alt_id",
    OBOINOWL + "consider": "consider",
    OBO_PURL + "IAO_0100001": "replaced_by",
}

SYNONYM_PRED2SCOPE: Dict[str, str] = {
    "hasExactSynonym": "EXACT",
    "hasBroadSynonym": "BROAD",
    "hasNarrowSynonym": "NARROW",
    "hasRelatedSynonym": "RELATED",
}


def get_obo_id(iri: str) -> str:
    """.../obo/GO_0008150 -> GO:0008150"""
    if iri.startswith(OBO_PURL):
        return iri[len(OBO_PURL):].replace("_", ":", 1)
    return iri


def get_fragment(iri: str) -> str:
    """.../obo/go#goslim_generic -> goslim_generic"""
    return iri.rsplit("#", 1)[-1]


def _obo_order(vals: List[str]) -> List[str]:
    """Sort as obo files do: case-insensitively."""
    return sorted(vals, key=lambda val: (val.lower(), val))


class OboGraphsReader(OBOReader):
    """OBOReader for OBO Graphs JSON: yields one GOTerm per term.

    Each node becomes the obo lines it stands for, read by OBOReader's record builder.
    """

    def __iter__(self) -> Iterator[GOTerm]:
        with open(self.obo_file) as fstream:
            graph: dict = json.load(fstream)["graphs"][0]
        self._init_header(graph.get("meta", {}))
        rel2tag: Dict[str, str] = self._init_typedefs(graph["nodes"])
        id2edges: Dict[str, List[str]] = self._get_id2edges(graph["edges"], rel2tag)
        terms: List[dict] = [n for n in graph["nodes"] if n.get("type") == "CLASS"]
        # Each alt_id is also a stub class; obo has it only as an alt_id line
        alt_ids: Set[str] = {
            pval["val"]
            for node in terms
            for pval in node.get("meta", {}).get("basicPropertyValues", [])
            if pval["pred"] == OBOINOWL + "hasAlternativeId"
        }
        for node in sorted(terms, key=lambda n: n["id"]):
            item_id: str = get_obo_id(node["id"])
            if item_id in alt_ids:
                continue
            rec = GOTerm(self.default_namespace)
            if self.optobj:
                self.optobj.init_datamembers(rec)
            for line in self._get_obo_lines(item_id, node, id2edges[item_id]):
                self._add_to_ref(rec, line)
            yield rec

    def _init_header(self, meta: dict) -> None:
        """Set format version, release, and default namespace."""
        pred2val: Dict[str, str] = {
            p["pred"]: p["val"] for p in meta.get("basicPropertyValues", [])
        }
        self.format_version = pred2val.get(OBOINOWL + "hasOBOFormatVersion")
        mtch = re.search(r"(releases/[^/]+)/", meta.get("version", ""))
        self.data_version = mtch.group(1) if mtch else None
        self.default_namespace = pred2val.get(
            OBOINOWL + "default-namespace", self.default_namespace
        )

    def _init_typedefs(self, nodes: List[dict]) -> Dict[str, str]:
        """Store relationship types as typedefs. Return relationship IRI -> obo name."""
        rel2tag: Dict[str, str] = {}
        for node in nodes:
            if node.get("type") != "PROPERTY":
                continue
            pred2val: Dict[str, str] = {
                p["pred"]: p["val"]
                for p in node.get("meta", {}).get("basicPropertyValues", [])
            }
            name = pred2val.get(OBOINOWL + "shorthand")
            if name is None:
                continue
            typedef = TypeDef()
            add_to_typedef(typedef, "id: " + name)
            add_to_typedef(typedef, "name: " + node.get("lbl", ""))
            add_to_typedef(
                typedef, "namespace: " + pred2val.get(OBOINOWL + "hasOBONamespace", "")
            )
            self.typedefs[name] = typedef
            rel2tag[node["id"]] = name
        return rel2tag

    @staticmethod
    def _get_id2edges(edges: List[dict], rel2tag: Dict[str, str]) -> Dict[str, List[str]]:
        """GO ID -> its is_a and relationship lines."""
        id2edges: Dict[str, List[str]] = cx.defaultdict(list)
        for edge in edges:
            pred, obj = edge["pred"], get_obo_id(edge["obj"])
            if pred == "is_a":
                id2edges[get_obo_id(edge["sub"])].append("is_a: " + obj)
            elif pred in rel2tag:
                line = "relationship: {R} {O}".format(R=rel2tag[pred], O=obj)
                id2edges[get_obo_id(edge["sub"])].append(line)
        return id2edges

    @staticmethod
    def _get_obo_lines(item_id: str, node: dict, edge_lines: List[str]) -> List[str]:
        """The obo lines this node stands for."""
        meta: dict = node.get("meta", {})
        lines: List[str] = ["id: " + item_id, "name: " + node.get("lbl", "")]
        for pval in meta.get("basicPropertyValues", []):
            tag = PRED2TAG.get(pval["pred"])
            if tag is not None:
                lines.append("{T}: {V}".format(T=tag, V=get_obo_id(pval["val"])))
        defn = meta.get("definition")
        if defn is not None:
            lines.append('def: "{V}" [{X}]'.format(
                V=defn.get("val", ""), X=", ".join(_obo_order(defn.get("xrefs", [])))))
        lines.extend("comment: " + comment for comment in meta.get("comments", []))
        lines.extend("subset: " + get_fragment(subset) for subset in meta.get("subsets", []))
        for syn in sorted(meta.get("synonyms", []), key=lambda s: (s["val"].lower(), s["val"])):
            typename = get_fragment(syn["synonymType"]) + " " if "synonymType" in syn else ""
            lines.append('synonym: "{V}" {S} {T}[{X}]'.format(
                V=syn["val"], S=SYNONYM_PRED2SCOPE[syn["pred"]], T=typename,
                X=", ".join(_obo_order(syn.get("xrefs", [])))))
        lines.extend("xref: " + xref["val"] for xref in meta.get("xrefs", []))
        lines.extend(edge_lines)
        if meta.get("deprecated"):
            lines.append("is_obsolete: true")
        return lines


# Copyright (C) 2010-present, H Tang et al. All rights reserved.
