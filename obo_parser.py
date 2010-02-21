#!/usr/bin/env python
# -*- coding: UTF-8 -*-

import os
import pprint
import sys
from itertools import chain
from exceptions import EOFError

typedef_tag, term_tag = "[Typedef]", "[Term]"


# stolen from http://kogs-www.informatik.uni-hamburg.de/~meine/python_tricks
def flatten(x):
    """flatten(sequence) -> list

    Returns a single, flat list which contains all elements retrieved
    from the sequence and all recursively contained sub-sequences
    (iterables).

    Examples:
    >>> [1, 2, [3,4], (5,6)]
    [1, 2, [3, 4], (5, 6)]
    >>> flatten([[[1,2,3], (42,None)], [4,5], [6], 7, MyVector(8,9,10)])
    [1, 2, 3, 42, None, 4, 5, 6, 7, 8, 9, 10]"""

    result = []
    for el in x:
        if hasattr(el, "__iter__") and not isinstance(el, basestring):
            result.extend(flatten(el))
        else:
            result.append(el)
    return result

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
            elif line.startswith("is_obsolete:") and line.split()[-1]=="true":
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
            if p.has_parent(term) or p.id==term:
                return True
        return False

    def get_all_parents(self):
        if not self.parents: return []
        return list(chain([[p.id] + p.get_all_parents()] \
                for p in self.parents))

    def has_child(self, term):
        for p in self.children:
            if p.has_child(term) or p.id==term:
                return True
        return False

    def get_all_children(self):
        if not self.children: return []
        return list(chain([[p.id] + p.get_all_children()] \
                for p in self.children))


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


def query_term(godag, term):
    rec = g[term]
    print >>sys.stderr, rec
    print >>sys.stderr, "all parents:", set(flatten(rec.get_all_parents()))
    print >>sys.stderr, "all children:", set(flatten(rec.get_all_children()))


def load_godag(obo_file="gene_ontology.1_2.obo"):

    g = GODag(obo_file)
    print >>sys.stderr, len(g), "nodes imported"

    return g

if __name__ == '__main__':

    import optparse
    p = optparse.OptionParser()
    p.add_option("-d", "--description", dest="desc", help="write term descriptions to stdout"
                 " from the obo file specified in args", action="store_true")
    p.add_option("-p", "--pickle", dest="make_pickle", 
                help="serialize the DAG data structure", action="store_true")

    (options, args) = p.parse_args()

    if not len(args):
        sys.exit(p.print_help())

    obo_file = args[0]
    assert os.path.exists(obo_file), "file %s not found!" % obo_file

    g = load_godag()

    if options.desc:
        g.write_dag()

    # run a test case
    query_term(g, "GO:0008135")


    # UNIMPLEMENTED
    """
    if options.make_pickle:
        import pickle
        pkl_file = obo_file + ".pkl"
        pickle.dump(g, file(pkl_file, "w"))
        print >>sys.stderr, "object serialized to", pkl_file
    """


