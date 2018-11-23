"""Get GO IDs from command-line arguments or from an ASCII file."""

from __future__ import print_function

__copyright__ = "Copyright (C) 2016-2019, DV Klopfenstein, H Tang. All rights reserved."
__author__ = "DV Klopfenstein"

import os
import re
from goatools.gosubdag.go_tasks import get_go2obj_unique
from goatools.godag.consts import Consts


class GetGOs(object):
    """Return a list of GO IDs for plotting."""

    def __init__(self, go2obj=None, max_gos=None):
        self.go2obj = go2obj
        self.max_gos = max_gos
        self.godagconsts = Consts()

    def get_goids(self, go_args, fin_goids, prt):
        """Return source GO IDs ."""
        goids = set()
        if fin_goids is not None:
            goids.update(self.rdtxt_gos(fin_goids, prt))
        if go_args:
            goids.update(self.get_goargs(go_args, prt))
        return goids

    def get_usrgos(self, fin_goids, prt):
        """Return source GO IDs ."""
        ret = self.get_goids(None, fin_goids, prt)
        # If there have been no GO IDs explicitly specified by the user
        if not ret:
            # If the GO-DAG is sufficiently small, print all GO IDs
            if self.max_gos is not None and len(self.go2obj) < self.max_gos:
                main_gos = set(o.id for go, o in self.go2obj.items() if go != o.id)
                go_leafs = set(go for go, o in self.go2obj.items() if not o.children)
                ret = go_leafs.difference(main_gos)
            else:
                raise RuntimeError("GO IDs NEEDED")
        go2obj = self.get_go2obj(ret)
        return get_go2obj_unique(go2obj)

    def get_go2obj(self, goids):
        """Return GO Terms for each user-specified GO ID. Note missing GO IDs."""
        goids = goids.intersection(self.go2obj.keys())
        if len(goids) != len(goids):
            goids_missing = goids.difference(goids)
            print("  {N} MISSING GO IDs: {GOs}".format(N=len(goids_missing), GOs=goids_missing))
        return {go:self.go2obj[go] for go in goids}

    @staticmethod
    def rdtxt_gos(go_file, prt):
        """Read GO IDs from a file."""
        goids_all = set()
        if not os.path.exists(go_file):
            raise RuntimeError("CAN NOT READ GO FILE: {FILE}\n".format(FILE=go_file))
        re_go = re.compile(r'(GO:\d{7})+?')
        re_com = re.compile(r'^\s*#')  # Lines starting with a '#' are comment lines and ignored
        with open(go_file) as ifstrm:
            for line in ifstrm:
                # Skip lines that are comments
                if re_com.search(line):
                    continue
                # Search for GO IDs on the line
                goids_found = re_go.findall(line)
                if goids_found:
                    goids_all.update(goids_found)
            if prt:
                prt.write("  {N} GO IDs READ: {TXT}\n".format(N=len(goids_all), TXT=go_file))
        return goids_all

    def get_goargs(self, go_args, prt):
        """Get GO IDs and colors for GO IDs from the GO ID runtime arguments."""
        goids = set()
        go2color = {}
        # Match on "GO ID" or "GO ID and color"
        re_gocolor = re.compile(r'(GO:\d{7})((?:#[0-9a-fA-F]{6})?)')
        for go_arg in go_args:
            mtch = re_gocolor.match(go_arg)
            if mtch:
                goid, color = mtch.groups()
                goids.add(goid)
                if color:
                    go2color[goid] = color
            elif go_arg in self.godagconsts.NS2GO:
                goids.add(self.godagconsts.NS2GO[go_arg])
            elif prt:
                prt.write("WARNING: UNRECOGNIZED ARG({})\n".format(go_arg))
        return goids

# Copyright (C) 2016-2019, DV Klopfenstein, H Tang. All rights reserved.
