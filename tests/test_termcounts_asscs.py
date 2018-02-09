"""Test TermCounts object used in Resnik and Lin similarity calculations."""

from __future__ import print_function

import os
import sys
from goatools.base import get_godag
from goatools.associations import dnld_assc
from goatools.semantic import TermCounts
from goatools.semantic import get_info_content

def test_semantic_similarity(usr_assc=None):
    """Computing basic semantic similarities between GO terms."""
    go2obj = get_go2obj()
    # goids = go2obj.keys()
    associations = [
        'gene_association.GeneDB_Lmajor',
        'gene_association.GeneDB_Pfalciparum',
        'gene_association.GeneDB_Tbrucei',
        'gene_association.GeneDB_tsetse',
        'gene_association.PAMGO_Atumefaciens',
        'gene_association.PAMGO_Ddadantii',
        #'gene_association.PAMGO_Mgrisea', # TBD Resolve DB_Name containing '|'
        'gene_association.PAMGO_Oomycetes',
        'gene_association.aspgd',
        'gene_association.cgd',
        'gene_association.dictyBase',
        'gene_association.ecocyc',
        'gene_association.fb',
        'gene_association.gonuts',
        #'gene_association.gramene_oryza', # DB_Name
        'gene_association.jcvi',
        'gene_association.mgi',
        'gene_association.pombase',
        'gene_association.pseudocap',
        'gene_association.reactome',
        'gene_association.rgd',
        'gene_association.sgd',
        'gene_association.sgn',
        'gene_association.tair',
        'gene_association.wb',
        'gene_association.zfin',
        'goa_chicken.gaf',
        'goa_chicken_complex.gaf',
        'goa_chicken_isoform.gaf',
        'goa_chicken_rna.gaf',
        'goa_cow.gaf',
        'goa_cow_complex.gaf',
        'goa_cow_isoform.gaf',
        'goa_cow_rna.gaf',
        'goa_dog.gaf',
        'goa_dog_complex.gaf',
        'goa_dog_isoform.gaf',
        'goa_dog_rna.gaf',
        'goa_human.gaf',
        'goa_human_complex.gaf',
        'goa_human_isoform.gaf',
        'goa_human_rna.gaf',
        'goa_pdb.gaf',
        'goa_pig.gaf',
        'goa_pig_complex.gaf',
        'goa_pig_isoform.gaf',
        'goa_pig_rna.gaf',
        'goa_uniprot_all.gaf',
        #'goa_uniprot_all_noiea.gaf',
    ]
    if usr_assc is not None:
        associations = [usr_assc]
    cwd = os.getcwd()
    for assc_name in associations:  # Limit test numbers for speed
        # Get all the annotations from arabidopsis.
        assc_gene2gos = dnld_assc(os.path.join(cwd, assc_name), go2obj, prt=None)

        # Calculate the information content of the single term, GO:0048364
        #       "Information content (GO:0048364) = 7.75481392334

        # First get the counts of each GO term.
        termcounts = TermCounts(go2obj, assc_gene2gos)
        go_cnt = termcounts.gocnts.most_common()
        #print termcounts.gocnts.most_common()

        if go_cnt:
            print("\n{ASSC}".format(ASSC=assc_name))
            print(sorted(termcounts.aspect_counts.most_common()))
            gocnt_max = go_cnt[0][1]
            prt_info(termcounts, go_cnt, None)
            prt_info(termcounts, go_cnt, gocnt_max/2.0)
            prt_info(termcounts, go_cnt, gocnt_max/10.0)

def prt_info(termcounts, go_cnt, max_val):
    """Print the information content of a frequently used GO ID."""
    go_id, cnt = get_goid(go_cnt, max_val)
    infocontent = get_info_content(go_id, termcounts)
    msg = 'Information content ({GO} {CNT:7,}) = {INFO:8.6f} {NAME}'
    print(msg.format(GO=go_id, CNT=cnt, INFO=infocontent, NAME=termcounts.go2obj[go_id].name))

def get_goid(go_cnt, max_val):
    """Get frequently used GO ID."""
    if max_val is not None:
        for goid, cnt in go_cnt:
            if cnt < max_val:
                return goid, cnt
        return go_cnt[-1][0], go_cnt[-1][1]
    return go_cnt[0][0], go_cnt[0][1]

def get_go2obj():
    """Read GODag and return go2obj."""
    godag = get_godag(os.path.join(os.getcwd(), "go-basic.obo"), loading_bar=None)
    return {go:o for go, o in godag.items() if not o.is_obsolete}

if __name__ == '__main__':
    ASSC_NAME = None if len(sys.argv) == 1 else sys.argv[1]
    test_semantic_similarity(ASSC_NAME)
