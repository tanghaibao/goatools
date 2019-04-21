#!/usr/bin/env python
"""Download evidenceontology/gaf-eco-mapping-derived.txt Save in eco2group.py"""

import os
# import sys
# from goatools.associations import dnld_ncbi_gene_file
# from goatools.anno.factory import get_objanno
from goatools.base import http_get

REPO = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..")


def dnld_eco2group():
    """Download evidenceontology/gaf-eco-mapping-derived.txt Save in eco2group.py"""
    fout_py = "goatools/anno/eco2group.py"
    url_ecomap = ("https://raw.githubusercontent.com/evidenceontology/evidenceontology/master/"
                  "gaf-eco-mapping-derived.txt")

    obj = _Run()
    fin_ecomap = obj.dnld(url_ecomap)
    eco2group = obj.read_ecomap(fin_ecomap)
    obj.wrpy(fout_py, eco2group)

# pylint: disable=superfluous-parens
class _Run(object):
    """Holds annotation filenames, downloads files if necessary"""

    REPO = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..")

    def wrpy(self, fout_py, eco2group):
        """Write eco2group to a Python module."""
        num_ecos = len(eco2group)
        with open(os.path.join(self.REPO, fout_py), 'w') as prt:
            prt.write('"""ECO ID to Evidence code group dictionary from evidenceontology"""\n\n')
            prt.write('# pylint: disable=too-many-lines\n')
            prt.write('# {N} ECO IDs\n'.format(N=num_ecos))
            prt.write('ECO2GRP = {\n')
            for eco, group in sorted(eco2group.items()):
                prt.write("    '{ECO}': '{GRP}',\n".format(ECO=eco, GRP=group))
            prt.write('}\n')
            print('{N} ECO IDs WROTE: {F}'.format(N=num_ecos, F=fout_py))

    @staticmethod
    def read_ecomap(fin_ecomap):
        """Read ECO mapping. Return dict, eco2group"""
        eco2group = {}
        ecos = set()
        with open(fin_ecomap) as ifstrm:
            for line in ifstrm:
                vals = line.split('\t')
                # print(vals)
                assert len(vals) >= 2
                eco = vals[0]
                grp = vals[1]
                len_grp = len(grp)
                assert len_grp == 2 or len_grp == 3, line
                ecos.add(eco)
                eco2group[eco] = grp
        num_keys = len(eco2group)
        print('{N} ECO IDs READ: {F}'.format(N=num_keys, F=fin_ecomap))
        assert num_keys == len(ecos)
        return eco2group

    def dnld(self, url_ecomap):
        """Get ECO mapping to group."""
        local_ecomap = os.path.join(self.REPO, url_ecomap.split('/')[-1])
        http_get(url_ecomap, local_ecomap)
        return local_ecomap


if __name__ == '__main__':
    dnld_eco2group()
