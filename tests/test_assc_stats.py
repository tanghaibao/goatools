"""Test association statistics."""

import sys
import os
from goatools.base import get_godag
from goatools.associations import dnld_assc
from goatools.associations import get_b2aset
from goatools.statsdescribe import StatsDescribe

def test_assc_stats(prt=sys.stdout):
    """Test association statistics."""
    associations = [
        ('hsa', 'goa_human.gaf'),         # human
        ('mus', 'gene_association.mgi'),  # mouse
        ('dme', 'gene_association.fb')]   # fly
    godag = get_godag(os.path.join(os.getcwd(), "go-basic.obo"), loading_bar=None)
    describe_go2obj(godag, prt)
    obj = StatsDescribe('Assc', "{:6,}")
    obj.prt_hdr(prt, "Assc.")
    for org, assc_name in associations:
        fin_assc = os.path.join(os.getcwd(), assc_name)
        describe_assc(org, fin_assc, godag, obj, prt)

def describe_go2obj(go2obj, prt):
    """Describe distribution of parent and child GO term counts."""
    # Related GO | # GO  | range    | 25th | median | 75th | mean | stddev
    # -----------|-------|----------|------|--------|------|------|-------
    # Parents    | 44961 | 0 to   8 |    1 |      1 |    2 |    2 |      1
    # Children   | 17597 | 1 to 480 |    1 |      2 |    4 |    4 |     10
    cnts_all = [(len(o.children), len(o.parents)) for go, o in go2obj.items() if go == o.id]
    cnts_c, cnts_p = zip(*cnts_all)
    cnts_c = [n for n in cnts_c if n != 0] # Remove leaf-level counts from reported stats
    cnts_p = [n for n in cnts_p if n != 0] # Remove top-level counts from reported stats
    obj = StatsDescribe('GO', "{:6,}")
    obj.prt_hdr(prt, "Related GO")
    obj.prt_data("Parents", cnts_p, prt)
    obj.prt_data("Children", cnts_c, prt)


def describe_assc(org, fin_assc, go2obj, obj, prt):
    """Report statistics for a single association."""
    # Assc.     | # Assc| range      | 25th | median | 75th | mean | stddev
    # ----------|-------|------------|------|--------|------|------|-------
    # hsa GOs   | 19394 | 1 to   212 |    5 |      9 |   17 |   13 |     14
    # hsa Genes | 17277 | 1 to 8,897 |    1 |      3 |    8 |   15 |    120
    #
    # mus GOs   | 19870 | 1 to   261 |    5 |     10 |   18 |   14 |     15
    # mus Genes | 17491 | 1 to 7,009 |    1 |      3 |    8 |   16 |    129
    #
    # dme GOs   | 12551 | 1 to   137 |    2 |      4 |    8 |    6 |      7
    # dme Genes |  7878 | 1 to 1,675 |    1 |      3 |    7 |   10 |     41
    gene2gos = dnld_assc(fin_assc, go2obj, prt=None) # Associations
    go2genes = get_b2aset(gene2gos)
    cnts_gos = [len(gos) for gos in gene2gos.values()]
    cnts_genes = [len(genes) for genes in go2genes.values()]
    obj.prt_data("{ORG} GOs".format(ORG=org), cnts_gos, prt)
    obj.prt_data("{ORG} Genes".format(ORG=org), cnts_genes, prt)

if __name__ == '__main__':
    test_assc_stats()
