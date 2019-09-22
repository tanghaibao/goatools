"""Test association statistics."""

import sys
import os
from goatools.base import get_godag
from goatools.associations import dnld_assc
from goatools.utils import get_b2aset
from goatools.statsdescribe import StatsDescribe


REPO = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..")

def test_assc_stats(prt=sys.stdout):
    """Test association statistics."""
    associations = [
        ('hsa', 'goa_human.gaf'), # human
        ('mus', 'mgi.gaf'),       # mouse
        ('dme', 'fb.gaf')]        # fly
    godag = get_godag(os.path.join(REPO, "go-basic.obo"), loading_bar=None)
    describe_go2obj(godag, prt)
    obj = StatsDescribe('Assc', "{:6,}")
    obj.prt_hdr(prt, "Assc.")
    for org, assc_name in associations:
        fin_assc = os.path.join(REPO, assc_name)
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
    # Assc.       | # Assc| range      | 25th | median | 75th | mean | stddev
    # ------------|-------|------------|------|--------|------|------|-------
    # hsa GO/gene | 19394 | 1 to   212 |    5 |      9 |   17 |   13 |     14
    # hsa gene/GO | 17277 | 1 to 8,897 |    1 |      3 |    8 |   15 |    120
    #
    # mus GO/gene | 19870 | 1 to   261 |    5 |     10 |   18 |   14 |     15
    # mus gene/GO | 17491 | 1 to 7,009 |    1 |      3 |    8 |   16 |    129
    #
    # dme GO/gene | 12551 | 1 to   137 |    2 |      4 |    8 |    6 |      7
    # dme gene/GO |  7878 | 1 to 1,675 |    1 |      3 |    7 |   10 |     41
    gene2gos = dnld_assc(fin_assc, go2obj, prt=None) # Associations
    go2genes = get_b2aset(gene2gos)
    assert gene2gos
    assert go2genes
    cnts_gos_p_gene = [len(gos) for gos in gene2gos.values()]
    cnts_genes_p_go = [len(genes) for genes in go2genes.values()]
    obj.prt_data("{ORG} GO/gene".format(ORG=org), cnts_gos_p_gene, prt)
    obj.prt_data("{ORG} gene/GO".format(ORG=org), cnts_genes_p_go, prt)

if __name__ == '__main__':
    test_assc_stats()
