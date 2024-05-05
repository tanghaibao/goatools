#!/usr/bin/env python
"""Tests writing GO hierarchy for BP, MF, CC."""
# https://github.com/tanghaibao/goatools/issues/163

import os
import sys

from goatools.base import get_godag
from goatools.associations import read_annotations
from goatools.semantic import TermCounts
from goatools.gosubdag.gosubdag import GoSubDag
from goatools.gosubdag.rpt.write_hierarchy import WrHierGO
from goatools.associations import dnld_ncbi_gene_file
from goatools.cli.wr_hierarchy import WrHierCli

REPO = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..")


def test_write_hier_bp_mf_cc():
    """Test that write hierarchy writes all: BP, MF, CC"""
    fin_anno = os.path.join(REPO, "gene2go")
    fin_dag = os.path.join(REPO, "go-basic.obo")
    _dnld_anno(fin_anno)

    print("\nTEST STORING ONLY ONE SPECIES")
    #### obj = Gene2GoReader(fin_anno)
    godag = get_godag(fin_dag)
    gene2gos = read_annotations(namespace="ALL")
    tcntobj = TermCounts(godag, gene2gos) if gene2gos else None
    gosubdag = GoSubDag(
        godag.keys(),
        godag,
        relationships=False,
        tcntobj=tcntobj,
        children=True,
        prt=sys.stdout,
    )

    print("Test using no_dup True: concise printing with no GO branches repeated")
    objwr = WrHierGO(gosubdag, no_dup=True)
    assert len(_wr_hier("nodup1", ["BP", "MF", "CC"], gosubdag.go2nt, objwr)) > 33000
    assert (
        len(
            _wr_hier(
                "nodup1",
                [
                    "BP",
                ],
                gosubdag.go2nt,
                objwr,
            )
        )
        > 25000
    )
    assert (
        len(
            _wr_hier(
                "nodup1",
                [
                    "MF",
                ],
                gosubdag.go2nt,
                objwr,
            )
        )
        > 10000
    )
    assert (
        len(
            _wr_hier(
                "nodup1",
                [
                    "CC",
                ],
                gosubdag.go2nt,
                objwr,
            )
        )
        > 4000
    )

    print("Test using no_dup False: verbose printing with duplicate copies branches")
    objwr = WrHierGO(gosubdag)
    # 2020 11:
    #     594,748 GO lines under GO:0008150
    #      23,199 GO lines under GO:0003674
    #       6,259 GO lines under GO:0005575
    #     624,206 items WROTE: tmp_test_wr_hier_BP_MF_CC.txt
    assert len(_wr_hier("nodup0", ["BP", "MF", "CC"], gosubdag.go2nt, objwr)) > 580000
    assert (
        len(
            _wr_hier(
                "nodup0",
                [
                    "BP",
                ],
                gosubdag.go2nt,
                objwr,
            )
        )
        > 500000
    )
    assert (
        len(
            _wr_hier(
                "nodup0",
                [
                    "MF",
                ],
                gosubdag.go2nt,
                objwr,
            )
        )
        > 20000
    )
    assert (
        len(
            _wr_hier(
                "nodup0",
                [
                    "CC",
                ],
                gosubdag.go2nt,
                objwr,
            )
        )
        > 5000
    )


def _wr_hier(desc, nss, go2nt, objwr):
    """Write hierarchy"""
    goids = WrHierCli.init_goids(nss, None, go2nt)
    fout_rpt = "tmp_test_wr_hier_{DESC}_{NSs}.txt".format(DESC=desc, NSs="_".join(nss))
    items_all = []
    with open(fout_rpt, "w", encoding="utf-8") as prt:
        for goid in goids:
            items_cur = objwr.prt_hier_down(goid, prt)
            items_all.extend(items_cur)
            print("{N:7,} GO lines under {GO}".format(N=len(items_cur), GO=goid))
    print("{N:7,} items WROTE: {TXT}".format(N=len(items_all), TXT=fout_rpt))
    return items_all


def _dnld_anno(file_anno):
    """Download the annotation file, if needed."""
    if os.path.exists(file_anno):
        assert os.path.getsize(file_anno) > 1000000, "BAD ANNO({F})".format(F=file_anno)
        return
    dnld_ncbi_gene_file(file_anno)
    assert os.path.isfile(file_anno), "MISSING ANNO({F})".format(F=file_anno)
    assert os.path.getsize(file_anno) > 1000000, "BAD ANNO({F})".format(F=file_anno)


if __name__ == "__main__":
    test_write_hier_bp_mf_cc()
