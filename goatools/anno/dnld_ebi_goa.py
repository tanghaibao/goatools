"""Download GOA files from the Gene Ontology Annotation (GOA) resource http://www.ebi.ac.uk/GOA."""

__copyright__ = "Copyright (C) 2016-2018, DV Klopfenstein, H Tang. All rights reserved."
__author__ = "DV Klopfenstein"

import os
import sys
from goatools.base import dnld_file

class DnldGoa(object):
    """Download files from the Gene Ontology Annotation (GOA) resource http://www.ebi.ac.uk/GOA."""

    # European Bioinformatics Institute (EMBL-EBI) ftp site
    ftp_pub = 'ftp://ftp.ebi.ac.uk/pub/'

    # Species available from ftp://ftp.ebi.ac.uk/pub/databases/GO/goa/
    species = [
        'arabidopsis',
        'chicken',
        'cow',
        'dicty',
        'dog',
        'fly',
        'human',
        'mouse',
        #'pdb',
        'pig',
        'rat',
        'uniprot',
        'worm',
        'yeast',
        'zebrafish',
    ]

    species_items = ['complex', 'isoform', 'rna']
    exts = ['gaf', 'gpa', 'gpi']

    def __init__(self):
        self.ftp_src_goa = os.path.join(self.ftp_pub, 'databases/GO/goa/')

    def dnld_goa(self, species, ext='gaf', item=None, fileout=None):
        """Download GOA source file name on EMBL-EBI ftp server."""
        basename = self.get_basename(species, ext, item)
        src = os.path.join(self.ftp_src_goa, species.upper(), "{F}.gz".format(F=basename))
        dst = os.path.join(os.getcwd(), basename) if fileout is None else fileout
        dnld_file(src, dst, prt=sys.stdout, loading_bar=None)
        return dst

    def get_basename(self, species, ext='gaf', item=None):
        """Get GOA basename for a specific species. Ex: goa_human.gaf"""
        assert ext in self.exts, " ".join(self.exts)
        if species == 'uniprot':
            species = 'uniprot_all' if item != 'gcrp' else 'uniprot_gcrp'
        if item is None:
            return 'goa_{SPECIES}.{EXT}'.format(SPECIES=species, EXT=ext)
        assert item in self.species_items
        return 'goa_{SPECIES}_{ITEM}.{EXT}'.format(SPECIES=species, ITEM=item, EXT=ext)


# Copyright (C) 2016-2018, DV Klopfenstein, H Tang. All rights reserved."
