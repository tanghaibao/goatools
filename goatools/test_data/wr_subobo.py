"""Write a small obo file to be used in tests.

   Create a small obo using a list of GO IDs as sources:

       import goatools.test_data.wr_subobo WrSubObo

       go_sources = ['GO:0003676', 'GO:0007516', 'GO:0036476']
       WrSubObo("go-basic.obo").wrobo("data/issue86.obo", go_sources)

"""

from __future__ import print_function

__copyright__ = "Copyright (C) 2010-2018, DV Klopfenstein, H Tang. All rights reserved."
__author__ = "DV Klopfenstein"

from goatools.obo_parser import GODag
from goatools.godag.go_tasks import CurNHigher


class WrSubObo(object):
    """Read a large GO-DAG from an obo file. Write a subset GO-DAG into a small obo file."""

    def __init__(self, fin_obo=None, optional_attrs=None, load_obsolete=None):
        self.fin_obo = fin_obo
        self.godag = GODag(fin_obo, optional_attrs, load_obsolete) if fin_obo is not None else None
        self.relationships = optional_attrs is not None and 'relationship' in optional_attrs

    def wrobo(self, fout_obo, goid_sources):
        """Write a subset obo file containing GO ID sources and their parents."""
        goids_all = self._get_goids_all(goid_sources)
        with open(fout_obo, 'w') as prt:
            self._prt_info(prt, goid_sources, goids_all)
            self.prt_goterms(prt, self.fin_obo, goids_all)
            print("  WROTE {N} GO TERMS: {OBO}\n".format(N=len(goids_all), OBO=fout_obo))

    @staticmethod
    def prt_goterms(fin_obo, goids, prt, b_prt=True):
        """Print the specified GO terms for GO IDs in arg."""
        b_trm = False
        with open(fin_obo) as ifstrm:
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
                        b_prt = line[4:14] in goids
                        if b_prt:
                            prt.write("[Term]\n")
                if b_prt:
                    prt.write(line)

    @staticmethod
    def get_goids(fin_obo, name):
        """Get GO IDs whose name matches given name."""
        goids = set()
        # pylint: disable=unsubscriptable-object
        goterm = None
        with open(fin_obo) as ifstrm:
            for line in ifstrm:
                if goterm is not None:
                    semi = line.find(':')
                    if semi != -1:
                        goterm[line[:semi]] = line[semi+2:].rstrip()
                    else:
                        if name in goterm['name']:
                            goids.add(goterm['id'])
                        goterm = None
                elif line[:6] == "[Term]":
                    goterm = {}
        return goids

    def _get_goids_all(self, go_sources):
        """Given GO ID sources and optionally the relationship attribute, return all GO IDs."""
        go2obj_user = {}
        objrel = CurNHigher(self.relationships, self.godag)
        objrel.get_id2obj_cur_n_high(go2obj_user, go_sources)
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
