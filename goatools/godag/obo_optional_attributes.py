"""Manage optional GO-DAG attributes."""

__copyright__ = "Copyright (C) 2015-2018, DV Klopfenstein, H Tang, All rights reserved."
__author__ = "DV Klopfenstein"

import re
import collections as cx


# pylint: disable=too-few-public-methods
class OboOptionalAttrs(object):
    """Manage optional GO-DAG attributes."""

    def __init__(self, optional_attrs=None):
        self.attr2cmp = self._init_compile_patterns(optional_attrs)
        self.optional_attrs = optional_attrs

    def update_rec(self, rec, line):
        """Update current GOTerm with optional record."""
        if 'def' in self.optional_attrs and line[:5] == "def: ":
            self._optattr_def(rec, line[5:])
        elif 'synonym' in self.optional_attrs and line[:9] == "synonym: ":
            self._optattr_synonym(rec, line[9:])
        elif 'relationship' in self.optional_attrs and line[:14] == "relationship: ":
            # relationships are now stored in a dict of sets. This mirrors
            # the structure implied by the GO DAG itself. The structure
            # that stores the relationships now looks likes this:
            #
            #  relationship = {
            #     'part_of': set(['GO:0021513', 'GO:0006310']),
            #     'negatively_regulates': set(['GO:0021910'])i
            # }
            self._optattr_relationship(rec, line[14:])
        elif 'xref' in self.optional_attrs and line[:6] == "xref: ":
            self._optattr_xref(rec, line[6:])
        elif 'subset' in self.optional_attrs and line[:8] == "subset: ":
            self._optattr_subset(rec, line[8:])
        elif 'comment' in self.optional_attrs and line[:9] == "comment: ":
            self._optattr_comment(rec, line[9:])

    @staticmethod
    def _optattr_def(rec, value):
        """Store optional attribute, 'defn' in a string."""
        # 'def' is a reserved word in python, do not use it as a Class attr.
        if not hasattr(rec, 'defn'):
            rec.defn = value
        else:
            raise Exception("ATTR(defn) ALREADY SET({VAL})".format(VAL=rec.defn))

    def _optattr_synonym(self, rec, value):
        """Store optional attribute, 'synonym' in a list."""
        if hasattr(rec, 'synonym'):
            rec.synonym.append(self._get_synonym(value))
        else:
            rec.synonym = [self._get_synonym(value)]

    def _optattr_relationship(self, rec, value):
        """Store optional attribute, 'relationship' in a dict of sets."""
        if hasattr(rec, 'relationship'):
            self._add_to_relationship(rec, value)
        else:
            rel, goid = value.split()[:2]
            rec.relationship = {rel: set([goid])}

    def _optattr_xref(self, rec, value):
        """Store optional attribute, 'xref' in a set."""
        if hasattr(rec, 'xref'):
            rec.xref.add(self._get_xref(value))
        else:
            rec.xref = set([self._get_xref(value)])

    @staticmethod
    def _optattr_subset(rec, value):
        """Store optional attribute, 'subset' in a set."""
        if hasattr(rec, 'subset'):
            rec.subset.add(value)
        else:
            rec.subset = set([value])

    @staticmethod
    def _optattr_comment(rec, value):
        """Store optional attribute, 'comment' in a string."""
        # 'def' is a reserved word in python, do not use it as a Class attr.
        if not hasattr(rec, 'comment'):
            rec.comment = value
        else:
            raise Exception("ATTR(comment) ALREADY SET({VAL})".format(VAL=rec.comment))

    def _get_synonym(self, line):
        """Given line, return optional attribute synonym value in a namedtuple."""
        # Example synonyms:
        # "peptidase inhibitor complex" EXACT [GOC:bf, GOC:pr]
        # "regulation of postsynaptic cytosolic calcium levels" EXACT syngo_official_label []
        # "tocopherol 13-hydroxylase activity" EXACT systematic_synonym []
        mtch = self.attr2cmp['synonym'].match(line)
        text, scope, typename, dbxrefs, _ = mtch.groups()
        typename = typename.strip()
        dbxrefs = set(dbxrefs.split(', ')) if dbxrefs else set()
        return self.attr2cmp['synonym nt']._make([text, scope, typename, dbxrefs])

    def _get_xref(self, line):
        """Given line, return optional attribute xref value in a dict of sets."""
        # Ex: xref      Wikipedia:Zygotene
        # Ex: Reactome:REACT_22295 "Addition of a third mannose to ..."
        mtch = self.attr2cmp['xref'].match(line)
        return mtch.group(1).replace(' ', '')

    @staticmethod
    def _add_to_relationship(rec, rel_value):
        """Add GO ID and its relationship to the optional 'relationship' data member."""
        rel, goid = rel_value.split()[:2]
        if rel in rec.relationship:
            rec.relationship[rel].add(goid)
        else:
            rec.relationship[rel] = set([goid])

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
