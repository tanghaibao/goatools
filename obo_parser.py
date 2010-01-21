#!/usr/bin/env python
# -*- coding: UTF-8 -*-

import os
from pprint import pprint
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
    __TODO__: this is not a parser, just a hack

    >>> reader = OBO_reader()
    >>> for rec in reader:
            print rec

    """

    def __init__(self, obo_file="gene_ontology.1_2.obo"):

        assert os.path.exists(obo_file), "%s not found!"%obo_file 
        self._handle = file(obo_file)
    
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
        self.id = "" 
        self.name = ""
        self.namespace = ""
        self.parents = set()

    def __str__(self):
        return "|".join([self.id, self.name, self.namespace])

    def __repr__(self):
        return self.__str__()

class GO_dag:

    def __init__(self, obo_file="gene_ontology.1_2.obo"):

        self.graph = {}
        self.load_obo_file(obo_file)

    def load_obo_file(self, obo_file):
        print "load obo file %s" % obo_file
        obo_reader = OBO_reader(obo_file)
        for rec in obo_reader:
            self.graph[rec.id] = rec
        #pprint(self.graph)

def write_desc(obo_file, out=sys.stdout):
    for rec in OBO_reader(obo_file):
        print >>out, rec.id + "\t" + rec.name 


if __name__ == '__main__':

    import optparse
    p = optparse.OptionParser()
    p.add_option("-d", dest="desc", help="write description file to stdout"
                 " from the obo file specified in args", action="store_true")

    opts, args = p.parse_args()

    if not (len(args) and opts.desc):
        p.print_help()
        my_dag = GO_dag()
        print len(my_dag.graph.keys()), "nodes imported"
        sys.exit()

    write_desc(args[0])

