"""Test that a Gene Ontology Enrichement Analysis can be run quietly"""

import os

from goatools.anno.idtogos_reader import IdToGosReader
from goatools.base import get_godag
from goatools.goea.go_enrichment_ns import GOEnrichmentStudyNS

__copyright__ = "Copyright (C) 2010-present, H Tang et al., All rights reserved."

REPO = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..")


def test_goea_quiet():
    """Test that a Gene Ontology Enrichement Analysis can be run quietly"""
    goeaobj = _get_goeaobj()
    study_fin = "{REPO}/tests/data/small_study".format(REPO=REPO)
    study_ids = [line.rstrip() for line in open(study_fin, encoding="utf-8")]
    print("\nTEST 1: GOEA run_study(study_ids)")
    goea_results1 = goeaobj.run_study(study_ids)
    print("{N} GOEA results for verbose GOEA".format(N=len(goea_results1)))

    print("\nTEST 2: GOEA run_study(study_ids, prt=None)")
    goea_results2 = goeaobj.run_study(study_ids, prt=None)
    print("{N} GOEA results for quiet GOEA".format(N=len(goea_results2)))

    # Original keyword is 'log'
    print("\nTEST 3: GOEA run_study(study_ids, log=None)")
    goea_results3 = goeaobj.run_study(study_ids, log=None)
    print("{N} GOEA results for quiet GOEA".format(N=len(goea_results3)))

    _chk_results(goea_results1, goea_results2)
    _chk_results(goea_results1, goea_results3)


def _get_goeaobj(methods=None):
    """Test GOEA with method, fdr."""
    # REad GODag
    obo_fin = os.path.join(REPO, "go-basic.obo")
    obo_dag = get_godag(obo_fin)
    # Read association
    fin_assc = "{REPO}/tests/data/small_association".format(REPO=REPO)
    objanno = IdToGosReader(fin_assc, godag=obo_dag)
    ns2assc = objanno.get_ns2assc()
    popul_fin = "{REPO}/tests/data/small_population".format(REPO=REPO)
    popul_ids = [line.rstrip() for line in open(popul_fin, encoding="utf-8")]
    goeaobj = GOEnrichmentStudyNS(popul_ids, ns2assc, obo_dag, methods=methods)
    return goeaobj


def _chk_results(results1, results2):
    """Check that results match"""
    # pylint: disable=line-too-long
    for res1, res2 in zip(results1, results2):
        assert res1.GO == res2.GO, "\nRES1: {R1}\nRES2: {R2}\n\n".format(
            R1=res1, R2=res2
        )
        assert (
            res1.p_bonferroni == res2.p_bonferroni
        ), "\nRES1: {R1}\nRES2: {R2}\n\n".format(R1=res1, R2=res2)
        assert res1.p_sidak == res2.p_sidak, "\nRES1: {R1}\nRES2: {R2}\n\n".format(
            R1=res1, R2=res2
        )
        assert res1.p_holm == res2.p_holm, "\nRES1: {R1}\nRES2: {R2}\n\n".format(
            R1=res1, R2=res2
        )


if __name__ == "__main__":
    test_goea_quiet()

# Copyright (C) 2010-present, H Tang et al., All rights reserved.
