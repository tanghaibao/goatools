"""Write a small obo file to be used in tests.

   Create a small obo using a list of GO IDs as sources:

       import goatools.test_data.wr_subobo WrSubObo

       go_sources = ['GO:0003676', 'GO:0007516', 'GO:0036476']
       WrSubObo("go-basic.obo").wrobo("data/issue86.obo", go_sources)

"""

__copyright__ = "Copyright (C) 2010-2018, DV Klopfenstein, H Tang. All rights reserved."
__author__ = "DV Klopfenstein"

import sys
from goatools.obo_parser import GODag
from goatools.godag.go_tasks import CurNHigher


# pylint: disable=too-few-public-methods
class WrSubObo(object):
    """Read a large GO-DAG from an obo file. Write a subset GO-DAG into a small obo file."""

    def __init__(self, fin_obo, optional_attrs=None, load_obsolete=None):
        self.fin_obo = fin_obo
        self.godag = GODag(fin_obo, optional_attrs, load_obsolete)
        self.relationships = optional_attrs is not None and 'relationship' in optional_attrs

    def wrobo(self, fout_obo, goid_sources):
        """Write a subset obo file containing GO ID sources and their parents."""
        goids_all = self._get_goids_all(goid_sources)
        b_trm = False
        b_prt = True
        with open(fout_obo, 'w') as prt:
            self._prt_info(prt, goid_sources, goids_all)
            with open(self.fin_obo) as ifstrm:
                for line in ifstrm:
                    if not b_trm:
                        if line[:6] == "[Term]":
                            b_trm = True
                            b_prt = False
                        elif line[:6] == "[Typedef]":
                            b_prt = True
                    else:
                        if line[:6] == 'id: GO':
                            b_trm = False
                            b_prt = line[4:14] in goids_all
                            if b_prt:
                                prt.write("[Term]\n")
                    if b_prt:
                        prt.write(line)
            sys.stdout.write("  WROTE {N} GO TERMS: {OBO}\n".format(N=len(goids_all), OBO=fout_obo))

    def _get_goids_all(self, go_sources):
        """Given GO ID sources and optionally the relationship attribute, return all GO IDs."""
        go2obj_user = {}
        objrel = CurNHigher(self.relationships, self.godag)
        objrel.get_go2obj_cur_n_high(go2obj_user, go_sources)
        goids = set(go2obj_user)
        for goterm in go2obj_user.values():
            if goterm.alt_ids:
                goids.update(goterm.alt_ids)
        return goids

    def _prt_info(self, prt, goid_sources, goids_all):
        """Print information describing how this obo setset was created."""
        prt.write("! Contains {N} GO IDs. Created using {M} GO sources:\n".format(
            N=len(goids_all), M=len(goid_sources)))
        for goid in goid_sources:
            prt.write("!    {GO}\n".format(GO=str(self.godag.get(goid, ""))))
        prt.write("\n")


# Copyright (C) 2010-2018, DV Klopfenstein, H Tang. All rights reserved.
