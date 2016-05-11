# Copyright 2010-2016 by Haibao Tang et al. All rights reserved.
#
# This code is part of the goatools distribution and goverend by its
# license. Please see the LICENSE file included with goatools.


"""Read and store Gene Ontology's obo file."""
# -*- coding: UTF-8 -*-
from __future__ import print_function
from collections import defaultdict
import sys
import os
import re

GraphEngines = ("pygraphviz", "pydot")

__copyright__ = "Copyright (C) 2010-2016, H Tang et al., All rights reserved."
__author__ = "various"

class OBOReader(object):
    """Read goatools.org's obo file. Load into this iterable class.

        Download obo from: http://geneontology.org/ontology/go-basic.obo

        >>> reader = OBOReader()
        >>> for rec in reader:
                print rec
    """

    def __init__(self, obo_file="go-basic.obo", optional_attrs=None):
        """Read obo file. Load dictionary."""
        self._init_optional_attrs(optional_attrs)
        self.format_version = None
        self.data_version = None
        self.typedefs = {}

        # True if obo file exists or if a link to an obo file exists.
        if os.path.isfile(obo_file):
            self.obo_file = obo_file
            # GOTerm attributes that are necessary for any operations:
        else:
            raise Exception("download obo file first\n "
                            "[http://geneontology.org/ontology/"
                            "go-basic.obo]")

    def __iter__(self):
        """Return one GO Term record at a time from an obo file."""
        # Written by DV Klopfenstein
        # Wait to open file until needed. Automatically close file when done.
        with open(self.obo_file) as fstream:
            rec_curr = None # Stores current GO Term
            typedef_curr = None  # Stores current typedef
            for lnum, line in enumerate(fstream):
                # obo lines start with any of: [Term], [Typedef], /^\S+:/, or /^\s*/
                if self.data_version is None:
                    self._init_obo_version(line)
                if line[0:6].lower() == "[term]":
                    rec_curr = self._init_goterm_ref(rec_curr, "Term", lnum)
                elif line[0:9].lower() == "[typedef]":
                    typedef_curr = self._init_typedef(rec_curr, "Typedef", lnum)
                elif rec_curr is not None or typedef_curr is not None:
                    line = line.rstrip() # chomp
                    if ":" in line:
                        if rec_curr is not None:
                            self._add_to_ref(rec_curr, line, lnum)
                        else:
                            self._add_to_typedef(typedef_curr, line, lnum)
                    elif line == "":
                        if rec_curr is not None:
                            yield rec_curr
                            rec_curr = None
                        elif typedef_curr is not None:
                            # Save typedef.
                            self.typedefs[typedef_curr.id] = typedef_curr
                            typedef_curr = None
                    else:
                        self._die("UNEXPECTED LINE CONTENT: {L}".format(L=line), lnum)
            # Return last record, if necessary
            if rec_curr is not None:
                yield rec_curr

    def _init_obo_version(self, line):
        """Save obo version and release."""
        if line[0:14] == "format-version":
            self.format_version = line[16:-1]
        if line[0:12] == "data-version":
            self.data_version = line[14:-1]

    def _init_goterm_ref(self, rec_curr, name, lnum):
        """Initialize new reference and perform checks."""
        if rec_curr is None:
            return GOTerm()
        msg = "PREVIOUS {REC} WAS NOT TERMINATED AS EXPECTED".format(REC=name)
        self._die(msg, lnum)

    def _init_typedef(self, typedef_curr, name, lnum):
        """Initialize new typedef and perform checks."""
        if typedef_curr is None:
            return TypeDef()
        msg = "PREVIOUS {REC} WAS NOT TERMINATED AS EXPECTED".format(REC=name)
        self._die(msg, lnum)

    def _add_to_ref(self, rec_curr, line, lnum):
        """Add new fields to the current reference."""
        # Written by DV Klopfenstein
        # Examples of record lines containing ':' include:
        #   id: GO:0000002
        #   name: mitochondrial genome maintenance
        #   namespace: biological_process
        #   def: "The maintenance of ...
        #   is_a: GO:0007005 ! mitochondrion organization
        mtch = re.match(r'^(\S+):\s*(\S.*)$', line)
        if mtch:
            field_name = mtch.group(1)
            field_value = mtch.group(2)
            if field_name == "id":
                self._chk_none(rec_curr.id, lnum)
                rec_curr.id = field_value
            elif field_name == "alt_id":
                rec_curr.alt_ids.append(field_value)
            elif field_name == "name":
                self._chk_none(rec_curr.name, lnum)
                rec_curr.name = field_value
            elif field_name == "namespace":
                self._chk_none(rec_curr.namespace, lnum)
                rec_curr.namespace = field_value
            elif field_name == "is_a":
                rec_curr._parents.append(field_value.split()[0])
            elif field_name == "is_obsolete" and field_value == "true":
                rec_curr.is_obsolete = True
            elif field_name in self.optional_attrs:
                self.update_rec(rec_curr, field_name, field_value)
        else:
            self._die("UNEXPECTED FIELD CONTENT: {L}\n".format(L=line), lnum)

    def update_rec(self, rec, name, value):
        """Update current GOTerm with optional record."""
        # 'def' is a reserved word in python, do not use it as a Class attr.
        if name == "def":
            name = "defn"

        # If we have a relationship, then we will split this into a further
        # dictionary.

        if hasattr(rec, name):
            if name not in self.attrs_scalar:
                if name not in self.attrs_nested:
                    getattr(rec, name).add(value)
                else:
                    self._add_nested(rec, name, value)
            else:
                raise Exception("ATTR({NAME}) ALREADY SET({VAL})".format(
                    NAME=name, VAL=getattr(rec, name)))
        else: # Initialize new GOTerm attr
            if name in self.attrs_scalar:
                setattr(rec, name, value)
            elif name not in self.attrs_nested:
                setattr(rec, name, set([value]))
            else:
                name = '_{:s}'.format(name)
                setattr(rec, name, defaultdict(list))
                self._add_nested(rec, name, value)

    def _add_to_typedef(self, typedef_curr, line, lnum):
        """Add new fields to the current typedef."""
        mtch = re.match(r'^(\S+):\s*(\S.*)$', line)
        if mtch:
            field_name = mtch.group(1)
            field_value = mtch.group(2).split('!')[0].rstrip()

            if field_name == "id":
                self._chk_none(typedef_curr.id, lnum)
                typedef_curr.id = field_value
            elif field_name == "name":
                self._chk_none(typedef_curr.name, lnum)
                typedef_curr.name = field_value
            elif field_name == "transitive_over":
                typedef_curr.transitive_over.append(field_value)
            elif field_name == "inverse_of":
                self._chk_none(typedef_curr.inverse_of, lnum)
                typedef_curr.inverse_of = field_value
            # Note: there are other tags that aren't imported here.
        else:
            self._die("UNEXPECTED FIELD CONTENT: {L}\n".format(L=line), lnum)

    def _add_nested(self, rec, name, value):
        """Adds a term's nested attributes."""
        # Remove comments and split term into typedef / target term.
        (typedef, target_term) = value.split('!')[0].rstrip().split(' ')

        # Save the nested term.
        getattr(rec, name)[typedef].append(target_term)

    def _init_optional_attrs(self, optional_attrs):
        """Prepare to store data from user-desired optional fields.

          Not loading these optional fields by default saves in space and speed.
          But allow the possibility for saving these fields, if the user desires,
            Including:
              comment consider def is_class_level is_metadata_tag is_transitive
              relationship replaced_by subset synonym transitive_over xref
        """
        # Written by DV Klopfenstein
        # Required attributes are always loaded. All others are optionally loaded.
        self.attrs_req = ['id', 'alt_id', 'name', 'namespace', 'is_a', 'is_obsolete']
        self.attrs_scalar = ['comment', 'defn',
                             'is_class_level', 'is_metadata_tag',
                             'is_transitive', 'transitive_over']
        self.attrs_nested = frozenset(['relationship'])
        # Allow user to specify either: 'def' or 'defn'
        #   'def' is an obo field name, but 'defn' is legal Python attribute name
        fnc = lambda aopt: aopt if aopt != "defn" else "def"
        if optional_attrs is None:
            optional_attrs = []
        elif isinstance(optional_attrs, str):
            optional_attrs = [fnc(optional_attrs)] if optional_attrs not in self.attrs_req else []
        elif isinstance(optional_attrs, list) or isinstance(optional_attrs, set):
            optional_attrs = set([fnc(f) for f in optional_attrs if f not in self.attrs_req])
        else:
            raise Exception("optional_attrs arg MUST BE A str, list, or set.")
        self.optional_attrs = optional_attrs


    def _die(self, msg, lnum):
        """Raise an Exception if file read is unexpected."""
        raise Exception("**FATAL {FILE}({LNUM}): {MSG}\n".format(
            FILE=self.obo_file, LNUM=lnum, MSG=msg))

    def _chk_none(self, init_val, lnum):
        """Expect these lines to be uninitialized."""
        if init_val is None or init_val is "":
            return
        self._die("FIELD IS ALREADY INITIALIZED", lnum)




class GOTerm:
    """
    GO term, actually contain a lot more properties than interfaced here
    """

    def __init__(self):
        self.id = ""                # GO:NNNNNNN
        self.name = ""              # description
        self.namespace = ""         # BP, CC, MF
        self._parents = []          # is_a basestring of parents
        self.parents = []           # parent records
        self.children = []          # children records
        self.level = None           # shortest distance from root node
        self.depth = None           # longest distance from root node
        self.is_obsolete = False    # is_obsolete
        self.alt_ids = []           # alternative identifiers

    def __str__(self):
        obsolete = "obsolete" if self.is_obsolete else ""
        return "{GO}\tlevel-{L:>02}\tdepth-{D:>02}\t{NAME} [{NS}] {OLD}".format(
            GO=self.id, L=self.level, D=self.depth,
            NAME=self.name, NS=self.namespace, OLD=obsolete)

    def __repr__(self):
        """Print GO id and all attributes in GOTerm class."""
        ret = ["GOTerm('{ID}'):".format(ID=self.id)]
        for key, val in self.__dict__.items():
            if isinstance(val, int) or isinstance(val, str):
                ret.append("{K}:{V}".format(K=key, V=val))
            else:
                ret.append("{K}: {V} items".format(K=key, V=len(val)))
                if len(val) < 10:
                    if not isinstance(val, dict):
                        for elem in val:
                            ret.append("  {ELEM}".format(ELEM=elem))
                    else:
                        for (typedef, terms) in val.items():
                            ret.append("  {TYPEDEF}: {NTERMS} items"
                                       .format(TYPEDEF=typedef,
                                               NTERMS=len(terms)))
                            for t in terms:
                                ret.append("    {TERM}".format(TERM=t))
        return "\n  ".join(ret)

    def has_parent(self, term):
        for p in self.parents:
            if p.id == term or p.has_parent(term):
                return True
        return False

    def has_child(self, term):
        for p in self.children:
            if p.id == term or p.has_child(term):
                return True
        return False

    def get_all_parents(self):
        all_parents = set()
        for p in self.parents:
            all_parents.add(p.id)
            all_parents |= p.get_all_parents()
        return all_parents

    def get_all_children(self):
        all_children = set()
        for p in self.children:
            all_children.add(p.id)
            all_children |= p.get_all_children()
        return all_children

    def get_all_parent_edges(self):
        all_parent_edges = set()
        for p in self.parents:
            all_parent_edges.add((self.id, p.id))
            all_parent_edges |= p.get_all_parent_edges()
        return all_parent_edges

    def get_all_child_edges(self):
        all_child_edges = set()
        for p in self.children:
            all_child_edges.add((p.id, self.id))
            all_child_edges |= p.get_all_child_edges()
        return all_child_edges

    def write_hier_rec(self, gos_printed, out=sys.stdout,
                      len_dash=1, max_depth=None, num_child=None, short_prt=False,
                      include_only=None, go_marks=None,
                      depth=1, dp="-"):
        """Write hierarchy for a GO Term record."""
        # Added by DV Klopfenstein
        GO_id = self.id
        # Shortens hierarchy report by only printing the hierarchy
        # for the sub-set of user-specified GO terms which are connected.
        if include_only is not None and GO_id not in include_only:
          return
        nrp = short_prt and GO_id in gos_printed
        if go_marks is not None:
          out.write('{} '.format('>' if GO_id in go_marks else ' '))
        if len_dash is not None:
            # Default character indicating hierarchy level is '-'.
            # '=' is used to indicate a hierarchical path printed in detail previously.
            letter = '-' if not nrp or not self.children else '='
            dp = ''.join([letter]*depth)
            out.write('{DASHES:{N}} '.format(DASHES=dp, N=len_dash))
        if num_child is not None:
            out.write('{N:>5} '.format(N=len(self.get_all_children())))
        out.write('{GO}\tL-{L:>02}\tD-{D:>02}\t{desc}\n'.format(
            GO=self.id, L=self.level, D=self.depth, desc=self.name))
        # Track GOs previously printed only if needed
        if short_prt:
          gos_printed.add(GO_id)
        # Do not print hierarchy below this turn if it has already been printed
        if nrp:
            return
        depth += 1
        if max_depth is not None and depth > max_depth:
            return
        for p in self.children:
            p.write_hier_rec(gos_printed, out, len_dash, max_depth, num_child, short_prt,
                include_only, go_marks,
                depth, dp)


class TypeDef(object):
    """
        TypeDef term. These contain more tags than included here, but these
        are the most important.
    """

    def __init__(self):
        self.id = ""                # GO:NNNNNNN
        self.name = ""              # description
        self.transitive_over = []   # List of other typedefs
        self.inverse_of = ""        # Name of inverse typedef.

    def __str__(self):
        ret = []
        ret.append("Typedef - {} ({}):".format(self.id, self.name))
        ret.append("  Inverse of: {}".format(self.inverse_of
                                             if self.inverse_of else "None"))
        if self.transitive_over:
            ret.append("  Transitive over:")
            for t in self.transitive_over:
                ret.append("    - {}".format(t))
        return "\n".join(ret)


class GODag(dict):

    def __init__(self, obo_file="go-basic.obo", optional_attrs=None):
        self.load_obo_file(obo_file, optional_attrs)

    def load_obo_file(self, obo_file, optional_attrs):

        print("load obo file %s" % obo_file, file=sys.stderr)
        reader = OBOReader(obo_file, optional_attrs)
        for rec in reader:
            self[rec.id] = rec
            for alt in rec.alt_ids:
                self[alt] = rec

        print("{OBO}: format-version({FMT}) data-version({REL})".format(
            OBO=obo_file, FMT=reader.format_version, REL=reader.data_version))

        # Save the typedefs and parsed optional_attrs
        self.typedefs = reader.typedefs
        self.optional_attrs = reader.optional_attrs

        self.populate_terms()
        print(len(self), "nodes imported", file=sys.stderr)

    def populate_terms(self):

        def _init_level(rec):
            if rec.level is None:
                if not rec.parents:
                    rec.level = 0
                else:
                    rec.level = min(_init_level(rec) for rec in rec.parents) + 1
            return rec.level

        def _init_depth(rec):
            if rec.depth is None:
                if not rec.parents:
                    rec.depth = 0
                else:
                    rec.depth = max(_init_depth(rec) for rec in rec.parents) + 1
            return rec.depth

        # Make parents and relationships references to the actual GO terms.
        for rec in self.values():
            rec.parents = [self[x] for x in rec._parents]

            if hasattr(rec, '_relationship'):
                rec.relationship = defaultdict(set)
                for (typedef, terms) in rec._relationship.items():
                    rec.relationship[typedef].update(set([self[x] for x in terms]))
                delattr(rec, '_relationship')

        # populate children, levels and add inverted relationships
        for rec in self.values():
            for p in rec.parents:
                if rec not in p.children:
                    p.children.append(rec)

            # Add invert relationships
            if hasattr(rec, 'relationship'):
                for (typedef, terms) in rec.relationship.items():
                    invert_typedef = self.typedefs[typedef].inverse_of
                    if invert_typedef:
                        # Add inverted relationship
                        for t in terms:
                            if not hasattr(t, 'relationship'):
                                t.relationship = defaultdict(set)
                            t.relationship[invert_typedef].add(rec)

            if rec.level is None:
                _init_level(rec)

            if rec.depth is None:
                _init_depth(rec)

    def write_dag(self, out=sys.stdout):
        """Write info for all GO Terms in obo file, sorted numerically."""
        for rec_id, rec in sorted(self.items()):
            print(rec, file=out)

    def write_hier_all(self, out=sys.stdout,
                      len_dash=1, max_depth=None, num_child=None, short_prt=False):
        """Write hierarchy for all GO Terms in obo file."""
        # Print: [biological_process, molecular_function, and cellular_component]
        for go_id in ['GO:0008150', 'GO:0003674', 'GO:0005575']:
          self.write_hier(go_id, out, len_dash, max_depth, num_child, short_prt, None)

    def write_hier(self, GO_id, out=sys.stdout,
                       len_dash=1, max_depth=None, num_child=None, short_prt=False,
                       include_only=None, go_marks=None):
        """Write hierarchy for a GO Term."""
        gos_printed = set()
        self[GO_id].write_hier_rec(gos_printed, out, len_dash, max_depth, num_child,
            short_prt, include_only, go_marks)

    @staticmethod
    def id2int(GO_id): return int(GO_id.replace("GO:", "", 1))

    def query_term(self, term, verbose=False):
        if term not in self:
            print("Term %s not found!" % term, file=sys.stderr)
            return

        rec = self[term]
        if verbose:
            print(rec)
            print("all parents:", rec.get_all_parents(), file=sys.stderr)
            print("all children:", rec.get_all_children(), file=sys.stderr)
        return rec

    def paths_to_top(self, term):
        """ Returns all possible paths to the root node

            Each path includes the term given. The order of the path is
            top -> bottom, i.e. it starts with the root and ends with the
            given term (inclusively).

            Parameters:
            -----------
            - term:
                the id of the GO term, where the paths begin (i.e. the
                accession 'GO:0003682')

            Returns:
            --------
            - a list of lists of GO Terms
        """
        # error handling consistent with original authors
        if term not in self:
            print("Term %s not found!" % term, file=sys.stderr)
            return

        def _paths_to_top_recursive(rec):
            if rec.level == 0:
                return [[rec]]
            paths = []
            for parent in rec.parents:
                top_paths = _paths_to_top_recursive(parent)
                for top_path in top_paths:
                    top_path.append(rec)
                    paths.append(top_path)
            return paths

        go_term = self[term]
        return _paths_to_top_recursive(go_term)

    def _label_wrap(self, label):
        wrapped_label = r"%s\n%s" % (label,
                                     self[label].name.replace(",", r"\n"))
        return wrapped_label

    def make_graph_pydot(self, recs, nodecolor,
                     edgecolor, dpi,
                     draw_parents=True, draw_children=True):
        """draw AMIGO style network, lineage containing one query record."""
        import pydot
        G = pydot.Dot(graph_type='digraph', dpi="{}".format(dpi)) # Directed Graph
        edgeset = set()
        usr_ids = [rec.id for rec in recs]
        for rec in recs:
            if draw_parents:
                edgeset.update(rec.get_all_parent_edges())
            if draw_children:
                edgeset.update(rec.get_all_child_edges())

        lw = self._label_wrap
        rec_id_set = set([rec_id for endpts in edgeset for rec_id in endpts])
        nodes = {str(ID):pydot.Node(
              lw(ID).replace("GO:",""),  # Node name
              shape="box",
              style="rounded, filled",
              # Highlight query terms in plum:
              fillcolor="beige" if ID not in usr_ids else "plum",
              color=nodecolor)
                for ID in rec_id_set}

        # add nodes explicitly via add_node
        for rec_id, node in nodes.items():
            G.add_node(node)

        for src, target in edgeset:
            # default layout in graphviz is top->bottom, so we invert
            # the direction and plot using dir="back"
            G.add_edge(pydot.Edge(nodes[target], nodes[src],
              shape="normal",
              color=edgecolor,
              label="is_a",
              dir="back"))

        return G

    def make_graph_pygraphviz(self, recs, nodecolor,
                     edgecolor, dpi,
                     draw_parents=True, draw_children=True):
        # draw AMIGO style network, lineage containing one query record
        import pygraphviz as pgv

        G = pgv.AGraph(name="GO tree")
        edgeset = set()
        for rec in recs:
            if draw_parents:
                edgeset.update(rec.get_all_parent_edges())
            if draw_children:
                edgeset.update(rec.get_all_child_edges())

        edgeset = [(self._label_wrap(a), self._label_wrap(b))
                   for (a, b) in edgeset]

        # add nodes explicitly via add_node
        # adding nodes implicitly via add_edge misses nodes
        # without at least one edge
        for rec in recs:
            G.add_node(self._label_wrap(rec.id))

        for src, target in edgeset:
            # default layout in graphviz is top->bottom, so we invert
            # the direction and plot using dir="back"
            G.add_edge(target, src)

        G.graph_attr.update(dpi="%d" % dpi)
        G.node_attr.update(shape="box", style="rounded,filled",
                           fillcolor="beige", color=nodecolor)
        G.edge_attr.update(shape="normal", color=edgecolor,
                           dir="back", label="is_a")
        # highlight the query terms
        for rec in recs:
            try:
                q = G.get_node(self._label_wrap(rec.id))
                q.attr.update(fillcolor="plum")
            except:
                continue

        return G

    def draw_lineage(self, recs, nodecolor="mediumseagreen",
                     edgecolor="lightslateblue", dpi=96,
                     lineage_img="GO_lineage.png", engine="pygraphviz",
                     gml=False, draw_parents=True, draw_children=True):
        assert engine in GraphEngines
        if engine == "pygraphviz":
            G = self.make_graph_pygraphviz(recs, nodecolor, edgecolor, dpi,
                              draw_parents=draw_parents, draw_children=draw_children)
        else:
            G = self.make_graph_pydot(recs, nodecolor, edgecolor, dpi,
                              draw_parents=draw_parents, draw_children=draw_children)

        if gml:
            import networkx as nx  # use networkx to do the conversion
            pf = lineage_img.rsplit(".", 1)[0]
            NG = nx.from_agraph(G) if engine == "pygraphviz" else nx.from_pydot(G)

            del NG.graph['node']
            del NG.graph['edge']
            gmlfile = pf + ".gml"
            nx.write_gml(NG, gmlfile)
            print("GML graph written to {0}".format(gmlfile), file=sys.stderr)

        print(("lineage info for terms %s written to %s" %
                             ([rec.id for rec in recs], lineage_img)), file=sys.stderr)

        if engine == "pygraphviz":
            G.draw(lineage_img, prog="dot")
        else:
            G.write_png(lineage_img)

    def update_association(self, association):
        bad_terms = set()
        for key, terms in list(association.items()):
            parents = set()
            for term in terms:
                try:
                    parents.update(self[term].get_all_parents())
                except:
                    bad_terms.add(term.strip())
            terms.update(parents)
        if bad_terms:
            print("terms not found: %s" % (bad_terms,), file=sys.stderr)

# Copyright (C) 2010-2016, H Tang et al., All rights reserved.
