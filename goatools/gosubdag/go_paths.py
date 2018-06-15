"""Collect GO paths."""

__copyright__ = "Copyright (C) 2016-2018, DV Klopfenstein, H Tang, All rights reserved."
__author__ = "DV Klopfenstein"

import collections as cx
import sys
from itertools import tee

class GoPaths(object):
    """Class for helping traverse GO paths."""

    adjdir = {
        None:  lambda go_obj: go_obj.parents + go_obj.children,
        True:  lambda go_obj: go_obj.parents,
        False: lambda go_obj: go_obj.children}

    # Adjacent direction can be any of:
    def get_paths_from_to(self, goobj_start, goid_end=None, dn0_up1=True):
        """Get a list of paths from goobj_start to either top or goid_end."""
        paths = []
        # Queue of terms to be examined (and storage for their paths)
        working_q = cx.deque([[goobj_start]])
        # Loop thru GO terms until we have examined all needed GO terms
        adjfnc = self.adjdir[dn0_up1]
        while working_q:
            #print "WORKING QUEUE LEN({})".format(len(working_q))
            path_curr = working_q.popleft()
            goobj_curr = path_curr[-1]
            go_adjlst = adjfnc(goobj_curr)
            #print 'END', goid_end, goobj_curr
            # If this GO term is the endpoint, Stop. Store path.
            if (goid_end is not None and goobj_curr.id == goid_end) or \
               (goid_end is None and not go_adjlst):
                paths.append(path_curr)
            # Else if this GO term is the not the end, add neighbors to path
            else:
                for go_neighbor in go_adjlst:
                    if go_neighbor not in path_curr:
                        #print "{}'s NEIGHBOR IS {}".format(goobj_curr.id, go_neighbor.id)
                        new_path = path_curr + [go_neighbor]
                        #sys.stdout.write(" {}'s {} {}\n".format(goobj_curr, up_dn, go_neighbor))
                        working_q.append(new_path)
        #self.prt_paths(paths)
        return paths

    @staticmethod
    def prt_paths(paths, prt=sys.stdout):
        """Print list of paths."""
        pat = "PATHES: {GO} L{L:02} D{D:02}\n"
        for path in paths:
            for go_obj in path:
                prt.write(pat.format(GO=go_obj.id, L=go_obj.level, D=go_obj.depth))
            prt.write("\n")

def get_paths_goobjs(go_objs, go_top=None, go2obj=None):
    """Given a list of GO objects, return: paths, user GOs as ints, all GO terms paths."""
    go_paths = []  # Collect all paths for go_objs
    go_all = set() # Collect all GO terms in all paths
    pathobj = GoPaths()
    for go_obj in go_objs:
        #print "?FIND PATHS FOR {}?".format(go_obj.id)
        if go_obj.id not in go_all: # GO not yet seen in paths already found
            #print "!FIND PATHS FOR {}!".format(go_obj.id)
            paths_curr = pathobj.get_paths_from_to(go_obj, go_top, True)
            if paths_curr:
                for path_goobjs in paths_curr:
                    for path_goobj in path_goobjs:
                        goid = path_goobj.id
                        if goid not in go_all:
                            go_all.add(goid)
                            go2obj[goid] = path_goobj
                # go_all.update(GO.id for path in paths_curr for GO in path)
                go_paths.extend(path for path in paths_curr)
    return go_paths, go_all

def paths2edges(paths):
    """[8079, 8135, 3723, 3676, 1901363, 5488, 3674] """
    edges_all = set()
    for path in paths:
        for edge in path2edges(path):
            edges_all.add(edge)
    return edges_all

def path2edges(path):
    """Given: [2000343, 32722, 1819]   Return: set([(2000343, 32722), (32722, 1819)])."""
    node_a, node_b = tee(path)
    next(node_b, None)
    return zip(node_a, node_b)

# Copyright (C) 2016-2018, DV Klopfenstein, H Tang, All rights reserved.
