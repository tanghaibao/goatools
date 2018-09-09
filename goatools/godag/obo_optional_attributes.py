"""Manage optional GO-DAG attributes."""

__copyright__ = "Copyright (C) 2015-2018, DV Klopfenstein, H Tang, All rights reserved."
__author__ = "DV Klopfenstein"

import re
import collections as cx


class OboOptionalAttrs(object):
    """Manage optional GO-DAG attributes."""

    attributes = set(['def', 'defn', 'synonym', 'relationship', 'xref', 'subset', 'comment'])

    def __init__(self, optional_attrs):
        assert optional_attrs
        self.optional_attrs = optional_attrs
        self.attr2cmp = self._init_compile_patterns(optional_attrs)

    def update_rec(self, rec, line):
        """Update current GOTerm with optional record."""
        if 'def' in self.optional_attrs and line[:5] == "def: ":
            assert not hasattr(rec, 'defn'), "ATTR(defn) ALREADY SET({VAL})".format(VAL=rec.defn)
            # Use 'defn' because 'def' is a reserved word in python
            rec.defn = line[5:]
        elif 'synonym' in self.optional_attrs and line[:9] == "synonym: ":
            rec.synonym.append(self._get_synonym(line[9:]))
        # http://geneontology.org/page/ontology-relations
        elif 'relationship' in self.optional_attrs and line[:14] == "relationship: ":
            # relationships are stored in a dict of sets, mirroring
            # the structure implied in the GO DAG. Example:
            #
            #  relationship = {
            #     'part_of': set(['GO:0021513', 'GO:0006310']),
            #     'regulates': set(['GO:0006313']),
            #     'negatively_regulates': set(['GO:0021910']),
            #     'positively_regulates': set(['GO:0006313']),
            # }
            rel, goid = line[14:].split()[:2]
            if rel not in rec.relationship:
                rec.relationship[rel] = set([goid])
            else:
                rec.relationship[rel].add(goid)
        elif 'xref' in self.optional_attrs and line[:6] == "xref: ":
            rec.xref.add(self._get_xref(line[6:]))
        elif 'subset' in self.optional_attrs and line[:8] == "subset: ":
            rec.subset.add(line[8:])
        elif 'comment' in self.optional_attrs and line[:9] == "comment: ":
            rec.comment = line[9:]

    def init_datamembers(self, rec):
        """Initialize current GOTerm with data members for storing optional attributes."""
        # pylint: disable=multiple-statements
        if 'synonym'      in self.optional_attrs: rec.synonym = []
        if 'xref'         in self.optional_attrs: rec.xref = set()
        if 'subset'       in self.optional_attrs: rec.subset = set()
        if 'comment'      in self.optional_attrs: rec.comment = ""
        if 'relationship' in self.optional_attrs:
            rec.relationship = {}
            rec.relationship_rev = {}

    def _get_synonym(self, line):
        """Given line, return optional attribute synonym value in a namedtuple.

        Example synonym and its storage in a namedtuple:
        synonym: "The other white meat" EXACT MARKETING_SLOGAN [MEAT:00324, BACONBASE:03021]
          text:     "The other white meat"
          scope:    EXACT
          typename: MARKETING_SLOGAN
          dbxrefs:  set(["MEAT:00324", "BACONBASE:03021"])

        Example synonyms:
          "peptidase inhibitor complex" EXACT [GOC:bf, GOC:pr]
          "regulation of postsynaptic cytosolic calcium levels" EXACT syngo_official_label []
          "tocopherol 13-hydroxylase activity" EXACT systematic_synonym []
        """
        mtch = self.attr2cmp['synonym'].match(line)
        text, scope, typename, dbxrefs, _ = mtch.groups()
        typename = typename.strip()
        dbxrefs = set(dbxrefs.split(', ')) if dbxrefs else set()
        return self.attr2cmp['synonym nt']._make([text, scope, typename, dbxrefs])

    def _get_xref(self, line):
        """Given line, return optional attribute xref value in a dict of sets."""
        # Ex: Wikipedia:Zygotene
        # Ex: Reactome:REACT_22295 "Addition of a third mannose to ..."
        mtch = self.attr2cmp['xref'].match(line)
        return mtch.group(1).replace(' ', '')

    @staticmethod
    def _init_compile_patterns(optional_attrs):
        """Compile search patterns for optional attributes if needed."""
        attr2cmp = {}
        if optional_attrs is None:
            return attr2cmp
        # "peptidase inhibitor complex" EXACT [GOC:bf, GOC:pr]
        # "blood vessel formation from pre-existing blood vessels" EXACT systematic_synonym []
        # "mitochondrial inheritance" EXACT []
        # "tricarboxylate transport protein" RELATED [] {comment="WIkipedia:Mitochondrial_carrier"}
        if 'synonym' in optional_attrs:
            attr2cmp['synonym'] = re.compile(r'"(\S.*\S)" ([A-Z]+) (.*)\[(.*)\](.*)$')
            attr2cmp['synonym nt'] = cx.namedtuple("synonym", "text scope typename dbxrefs")
        # Wikipedia:Zygotene
        # Reactome:REACT_27267 "DHAP from Ery4P and PEP, Mycobacterium tuberculosis"
        if 'xref' in optional_attrs:
            attr2cmp['xref'] = re.compile(r'^(\S+:\s*\S+)\b(.*)$')
        return attr2cmp

    @staticmethod
    def get_optional_attrs(optional_attrs):
        """Prepare to store data from user-desired optional fields.

          Not loading these optional fields by default saves in space and speed.
          But allow the possibility for saving these fields, if the user desires,
            Including:
              comment consider def is_class_level is_metadata_tag is_transitive
              relationship replaced_by subset synonym transitive_over xref
        """
        attrs_opt = set(['def', 'defn', 'synonym', 'relationship', 'xref', 'subset', 'comment'])
        # Required attributes are always loaded. All others are optionally loaded.
        # Allow user to specify either: 'def' or 'defn'
        #   'def' is an obo field name, but 'defn' is legal Python attribute name
        getnm = lambda aopt: aopt if aopt != "defn" else "def"
        # pylint: disable=redefined-variable-type
        opts = None
        if isinstance(optional_attrs, str) and optional_attrs in attrs_opt:
            opts = set([getnm(optional_attrs)])
        else:
            opts = set([getnm(f) for f in optional_attrs if f in attrs_opt])
        if opts:
            return opts


# Copyright (C) 2015-2018, DV Klopfenstein, H Tang, All rights reserved.
