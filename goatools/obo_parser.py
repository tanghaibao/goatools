#!/usr/bin/env python
# -*- coding: UTF-8 -*-

import os
import pprint
import sys
from itertools import chain
from exceptions import EOFError

typedef_tag, term_tag = "[Typedef]", "[Term]"

def after_colon(line):
    # macro for getting anything after the :
    return line.split(":", 1)[1].strip()

def read_until(handle, start):
    # read each line until it has a certain start, and then puts the start tag back
    while 1:
        pos = handle.tell()
        line = handle.readline()
        if not line: 
            break
        if line.startswith(start):
            handle.seek(pos)
            return
    raise EOFError, "%s tag cannot be found"


class OBOReader:
    """
    parse obo file, usually the most updated can be downloaded from
    http://www.geneontology.org/ontology/obo_format_1_2/gene_ontology.1_2.obo

    >>> reader = OBOReader()
    >>> for rec in reader:
            print rec

    """

    def __init__(self, obo_file="gene_ontology.1_2.obo"):

        try:
            self._handle = file(obo_file)
        except:
            print >>sys.stderr, \
                "download obo file first\n " \
                "[http://geneontology.org/ontology/obo_format_1_2/gene_ontology.1_2.obo]"
            sys.exit(1)
    
    def __iter__(self):

        term_tag = "[Term]"
        line = self._handle.readline()
        if not line.startswith(term_tag):
            read_until(self._handle, term_tag)
        while 1:
            yield self.next()

    def next(self):

        lines = []
        line = self._handle.readline()
        if not line or line.startswith(typedef_tag):
            raise StopIteration

        # read until the next tag and save everything in between
        while 1:
            pos = self._handle.tell() # save current postion for roll-back
            line = self._handle.readline()
            if line.startswith(typedef_tag) or line.startswith(term_tag):
                self._handle.seek(pos) # roll-back
                break
            lines.append(line)

        rec = GOTerm()
        for line in lines:
            if line.startswith("id:"):
                rec.id = after_colon(line) 
            elif line.startswith("name:"):
                rec.name = after_colon(line) 
            elif line.startswith("namespace:"):
                rec.namespace = after_colon(line) 
            elif line.startswith("is_a:"):
                rec._parents.append(after_colon(line).split()[0])
            elif line.startswith("is_obsolete:") and after_colon(line)=="true":
                rec.is_obsolete = True

        return rec


class GOTerm:
    """
    GO term, actually contain a lot more properties than interfaced here
    """

    def __init__(self):
        self.id = ""             # GO:xxxxxx
        self.name = ""           # description
        self.namespace = ""      # BP, CC, MF
        self._parents = []       # is_a basestring of parents
        self.parents  = []       # parent records
        self.children = []       # children records
        self.level = -1          # distance from root node
        self.is_obsolete = False # is_obsolete

    def __str__(self):
        obsolete = "obsolete" if self.is_obsolete else ""
        return "%s\tlevel-%02d\t%s [%s] %s" % \
                    (self.id, self.level, self.name,
                     self.namespace, obsolete) 

    def __repr__(self):
        return "GOTerm('%s')" % (self.id) 

    def has_parent(self, term):
        for p in self.parents:
            if p.id==term or p.has_parent(term):
                return True
        return False

    def has_child(self, term):
        for p in self.children:
            if p.id==term or p.has_child(term):
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


class GODag:

    def __init__(self, obo_file="gene_ontology.1_2.obo"):

        self.graph = {}
        self.load_obo_file(obo_file)

    def __len__(self):
        return len(self.graph)

    def __getitem__(self, key):
        return self.graph[key]

    def keys(self):
        return self.graph.keys()

    def load_obo_file(self, obo_file):
        
        print >>sys.stderr, "load obo file %s" % obo_file
        obo_reader = OBOReader(obo_file)
        for rec in obo_reader:
            self.graph[rec.id] = rec

        self.populate_terms()
        print >>sys.stderr, len(self), "nodes imported"

    def populate_terms(self):
        
        def depth(rec):
            if rec.level < 0:
                if not rec.parents:
                    rec.level = 0
                else:
                    rec.level = min(depth(rec) for rec in rec.parents) + 1
            return rec.level

        # make the parents references to the GO terms
        for rec in self.graph.itervalues():
            rec.parents = [self[x] for x in rec._parents]

        # populate children and levels
        for rec in self.graph.itervalues():
            for p in rec.parents:
                p.children.append(rec)

            if rec.level < 0: 
                depth(rec)


    def write_dag(self, out=sys.stdout):

        for rec_id, rec in sorted(self.graph.items()):
            print >>out, rec


    def query_term(self, term, draw_lineage=False):
        try:
            rec = self[term]
        except:
            print >>sys.stderr, "Term %s not found!" % term
            return
        print >>sys.stderr, rec
        print >>sys.stderr, "all parents:", rec.get_all_parents()
        print >>sys.stderr, "all children:", rec.get_all_children()
        if draw_lineage:
            self.draw_lineage(rec)


    def draw_lineage(self, rec):
        # draw AMIGO style network, lineage containing one query record
        try:
            import pydot
        except:
            print >>sys.stderr, "pydot not installed, lineage not drawn!"
            print >>sys.stderr, "try `easy_install pydot`"
            return
        
        G = pydot.Dot() 
        edgeset = rec.get_all_parent_edges() | rec.get_all_child_edges()
        for src, target in edgeset:
            # ":" is interpreted as something else in pydot
            src, target = src.replace(":","_"), target.replace(":", "_")
            G.add_edge(pydot.Edge(src, target))

        lineage_img = "%s.jpg" % rec.id.replace(":", "_")
        print >>sys.stderr, "lineage info for term %s written to %s" %\
                (rec.id, lineage_img)


if __name__ == '__main__':

    import optparse
    p = optparse.OptionParser()
    p.add_option("-d", "--description", dest="desc", 
            help="write term descriptions to stdout" \
                 " from the obo file specified in args", action="store_true")
    p.add_option("-t", "--term", dest="term", help="write the parents and children" \
            "of the query term", action="store", type="string", default=None)

    (options, args) = p.parse_args()

    if not len(args):
        sys.exit(p.print_help())

    obo_file = args[0]
    assert os.path.exists(obo_file), "file %s not found!" % obo_file

    g = GODag(obo_file)

    if options.desc:
        g.write_dag()

    # run a test case
    if options.term is not None:
        g.query_term(options.term, draw_lineage=True)

