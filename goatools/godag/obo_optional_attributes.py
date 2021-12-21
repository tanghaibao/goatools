"""Manage optional GO-DAG attributes."""

__copyright__ = (
    "Copyright (C) 2015-present, DV Klopfenstein, H Tang, All rights reserved."
)
__author__ = "DV Klopfenstein"

import re
import collections as cx


class OboOptionalAttrs:
    """Manage optional GO-DAG attributes."""

    optional_exp = {
        "def",
        "defn",
        "synonym",
        "relationship",
        "xref",
        "subset",
        "comment",
        "consider",
        "replaced_by",
    }

    def __init__(self, optional_attrs):
        assert optional_attrs
        self.optional_attrs = optional_attrs.intersection(self.optional_exp)
        self.fncs_inirec = self._init_fncs_inirec()
        self.fncs = self._init_fncs()

    def update_rec(self, rec, line):
        """Update current GOTerm with optional record."""
        for fnc_chkline, fnc_updaterec in self.fncs:
            if fnc_chkline(line):
                fnc_updaterec(rec, line)

    def _init_fncs(self):
        """Initialize functions to check for optional attributes and update GOTerm"""
        fncs = set()
        optional_attrs = self.optional_attrs
        if "def" in optional_attrs:
            fncs.add(self._get_fncs_def())
        if "synonym" in optional_attrs:
            fncs.add(self._get_fncs_synonym())
        if "relationship" in optional_attrs:
            fncs.add(self._get_fncs_relationship())
        if "xref" in optional_attrs:
            fncs.add(self._get_fncs_xref())
        if "subset" in optional_attrs:
            fncs.add(self._get_fncs_subset())
        if "comment" in optional_attrs:
            fncs.add(self._get_fncs_comment())
        if "consider" in optional_attrs:
            fncs.add(self._get_fncs_consider())
        if "replaced_by" in optional_attrs:
            fncs.add(self._get_fncs_replaced_by())
        return fncs

    @staticmethod
    def _get_fncs_def():
        def fnc_chkline(line):
            return line[:5] == "def: "

        def fnc_updaterec(rec, line):
            assert not hasattr(rec, "defn"), "ATTR(defn) ALREADY SET({VAL})".format(
                VAL=rec.defn
            )
            # Use 'defn' because 'def' is a reserved word in python
            rec.defn = line[5:]

        return fnc_chkline, fnc_updaterec

    @staticmethod
    def _get_fncs_synonym():
        def fnc_chkline(line):
            return line[:9] == "synonym: "

        def fnc_updaterec(rec, line):
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
            mtch = fnc_updaterec.cmpd.match(line[9:])
            text, scope, typename, dbxrefs, _ = mtch.groups()
            typename = typename.strip()
            dbxrefs = set(dbxrefs.split(", ")) if dbxrefs else set()
            ntd = fnc_updaterec.ntobj._make([text, scope, typename, dbxrefs])
            rec.synonym.append(ntd)

        fnc_updaterec.cmpd = re.compile(r'"(\S.*\S)" ([A-Z]+) (.*)\[(.*)\](.*)$')
        fnc_updaterec.ntobj = cx.namedtuple("synonym", "text scope typename dbxrefs")
        return fnc_chkline, fnc_updaterec

    @staticmethod
    def _get_fncs_relationship():
        # http://geneontology.org/page/ontology-relations
        def fnc_chkline(line):
            return line[:14] == "relationship: "

        def fnc_updaterec(rec, line):
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
                rec.relationship[rel] = {goid}
            else:
                rec.relationship[rel].add(goid)

        return fnc_chkline, fnc_updaterec

    @staticmethod
    def _get_fncs_xref():
        def fnc_chkline(line):
            return line[:6] == "xref: "

        def fnc_updaterec(rec, line):
            s = line[6:]
            match = fnc_updaterec.cmpd.match(s)
            if match:
                s = match.group(1)
                rec.xref.add(s.replace(" ", ""))
            else:  # occasionally just 'EC 2.7.1.190'
                rec.xref.add(s.replace(" ", ":"))

        fnc_updaterec.cmpd = re.compile(r"^(\S+:\s*\S+)\b(.*)$")
        return fnc_chkline, fnc_updaterec

    @staticmethod
    def _get_fncs_subset():
        def fnc_chkline(line):
            return line[:8] == "subset: "

        def fnc_updaterec(rec, line):
            rec.subset.add(line[8:])

        return fnc_chkline, fnc_updaterec

    @staticmethod
    def _get_fncs_comment():
        def fnc_chkline(line):
            return line[:9] == "comment: "

        def fnc_updaterec(rec, line):
            rec.comment = line[9:]

        return fnc_chkline, fnc_updaterec

    @staticmethod
    def _get_fncs_consider():
        """Get optional attribute functions"""

        def fnc_chkline(line):
            return line[:10] == "consider: "

        def fnc_updaterec(rec, line):
            rec.consider.add(line[10:])

        return fnc_chkline, fnc_updaterec

    @staticmethod
    def _get_fncs_replaced_by():
        def fnc_chkline(line):
            return line[:13] == "replaced_by: "

        def fnc_updaterec(rec, line):
            rec.replaced_by = line[13:]

        return fnc_chkline, fnc_updaterec

    def init_datamembers(self, rec):
        """Initialize current GOTerm with data members for storing optional attributes."""
        for fnc_ini in self.fncs_inirec:
            fnc_ini(rec)

    @staticmethod
    def _init_synonym(rec):
        rec.synonym = []

    @staticmethod
    def _init_xref(rec):
        rec.xref = set()

    @staticmethod
    def _init_subset(rec):
        rec.subset = set()

    @staticmethod
    def _init_comment(rec):
        rec.comment = ""

    @staticmethod
    def _init_replaced_by(rec):
        rec.replaced_by = ""

    @staticmethod
    def _init_consider(rec):
        rec.consider = set()

    @staticmethod
    def _init_relationship(rec):
        rec.relationship = {}
        rec.relationship_rev = {}

    def _init_fncs_inirec(self):
        """Initialize functions to check for optional attributes and update GOTerm"""
        fncs = set()
        optional_attrs = self.optional_attrs
        if "synonym" in optional_attrs:
            fncs.add(self._init_synonym)
        if "relationship" in optional_attrs:
            fncs.add(self._init_relationship)
        if "xref" in optional_attrs:
            fncs.add(self._init_xref)
        if "subset" in optional_attrs:
            fncs.add(self._init_subset)
        if "comment" in optional_attrs:
            fncs.add(self._init_comment)
        if "consider" in optional_attrs:
            fncs.add(self._init_consider)
        if "replaced_by" in optional_attrs:
            fncs.add(self._init_replaced_by)
        return fncs

    @staticmethod
    def get_optional_attrs(optional_attrs, attrs_opt):
        """Prepare to store data from user-desired optional fields.

        Not loading these optional fields by default saves in space and speed.
        But allow the possibility for saving these fields, if the user desires,
          Including:
            comment consider def is_class_level is_metadata_tag is_transitive
            relationship replaced_by subset synonym transitive_over xref
        """
        # Required attributes are always loaded. All others are optionally loaded.
        # Allow user to specify either: 'def' or 'defn'
        #   'def' is an obo field name, but 'defn' is legal Python attribute name
        try:
            iter(optional_attrs)
        except TypeError:
            pat = (
                "**FATAL: GODag's optional_attrs MUST BE A SET CONTAINING ANY OF: {ATTRS}\n"
                "           "
                "**FATAL: BAD GODag optional_attrs({BADVAL})"
            )
            msg = pat.format(ATTRS=" ".join(attrs_opt), BADVAL=optional_attrs)
            raise TypeError(msg)
        getnm = lambda aopt: aopt if aopt != "defn" else "def"
        opts = None
        if isinstance(optional_attrs, str) and optional_attrs in attrs_opt:
            opts = {getnm(optional_attrs)}
        else:
            opts = set(getnm(f) for f in optional_attrs if f in attrs_opt)
        return opts


# Copyright (C) 2015-present, DV Klopfenstein, H Tang, All rights reserved.
