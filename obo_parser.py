#!/usr/bin/env python
# -*- coding: UTF-8 -*-

import os
import pprint
import sys
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


class OBO_reader:
    """
    parse obo file, usually the most updated can be downloaded from
    http://www.geneontology.org/ontology/obo_format_1_2/gene_ontology.1_2.obo

    >>> reader = OBO_reader()
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

        rec = GO_term()
        for line in lines:
            if line.startswith("id:"):
                rec.id = after_colon(line) 
            elif line.startswith("name:"):
                rec.name = after_colon(line) 
            elif line.startswith("namespace:"):
                rec.namespace = after_colon(line) 
            elif line.startswith("is_a:"):
                rec.parents.add(after_colon(line).split()[0])

        return rec


class GO_term:
    """
    GO term, actually contain a lot more properties than interfaced here
    """

    def __init__(self):
        self.id = ""            # GO:xxxxxx
        self.name = ""          # description
        self.namespace = ""     # BP, CC, MF
        self.parents = set()    # is_a
        self.level = -1         # distance from root node

    def __str__(self):

        return "%s\tlevel-%02d\t%s [%s]" % \
                    (self.id, self.level, self.name, self.namespace) 

    __repr__ = __str__

    def getlevel(self):
        # TODO: I think what I am doing is reasonable here
        # taking the smallest level it can be... need to check though
        # the confusing part is that the term can be at different levels
        # depending on what ontology lineage you are tracking
        # see for example GO:0000022 [search on AMIGO and look at the graph]
        # it can be 4 or 7, and I am calling it 4

        if self.level < 0: # uninitialized 

            if not self.parents: # root or obsolete terms
                self.level = 0
            else:
                self.level = min(p.getlevel() \
                        for p in self.parents) + 1

        return self.level


class GO_dag:

    def __init__(self, obo_file="gene_ontology.1_2.obo"):

        self.graph = {}
        self.load_obo_file(obo_file)

    def __len__(self):
        return len(self.graph)

    def __getitem__(self, key):

        return self.graph[key]

    def load_obo_file(self, obo_file):
        
        print >>sys.stderr, "load obo file %s" % obo_file
        obo_reader = OBO_reader(obo_file)
        for rec in obo_reader:
            self.graph[rec.id] = rec

        self.populate_terms()

    def populate_terms(self):

        # make the parents references to the GO terms
        for rec in self.graph.values():
            rec.parents = [self[key] for key in rec.parents]

        # populate levels
        for rec in self.graph.values():
            if rec.level < 0:
                rec.getlevel()

    def write_dag(self, out=sys.stdout):

        for rec_id, rec in sorted(self.graph.items()):
            print >>out, rec


if __name__ == '__main__':

    import optparse
    p = optparse.OptionParser()
    p.add_option("-d", dest="desc", help="write description file to stdout"
                 " from the obo file specified in args", action="store_true")

    options, args = p.parse_args()

    if not len(args):
        sys.exit(p.print_help())

    obo_file = args[0]
    assert os.path.exists(obo_file), "file %s not found!" % obo_file

    g = GO_dag(obo_file)
    print >>sys.stderr, len(g), "nodes imported"

    if options.desc:
        g.write_dag()


