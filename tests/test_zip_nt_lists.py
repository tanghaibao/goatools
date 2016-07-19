"""Test combining lists whose elements are namedtuples or class objects."""

import os
import sys
import collections as cx
from goatools.wr_tbl import zip_nt_lists
from goatools.obo_parser import GODag
from goatools.go_enrichment import GOEnrichmentStudy, get_goea_nts_all, get_goea_nts
from goatools.associations import read_associations

def test_zip_nt_lists(prt=sys.stdout):
    """Test combining lists whose elements are namedtuples or class objects."""
    ntobj = cx.namedtuple("Nt", "idx")
    obo_dag, goea_results = get_goea_results()
    # Zip a list of namedtuples and another list of namedtuples
    goea_nts = get_goea_nts(goea_results)
    lst2_nts = [ntobj._make([i]) for i in range(len(goea_nts))]
    # Combine lists into a single list whose elements are a namedtuple
    lst_all = zip_nt_lists([lst2_nts, goea_nts])
    assert lst_all[0]._fields == lst2_nts[0]._fields + goea_nts[0]._fields

def get_goea_results(method="fdr_bh", log=sys.stdout):
    """Get GOEA results."""
    ROOT = os.path.dirname(os.path.abspath(__file__)) + "/data/"
    obo_dag = GODag(ROOT + "goslim_generic.obo")
    assoc = read_associations(ROOT + "slim_association", no_top=True)
    popul_ids = [line.rstrip() for line in open(ROOT + "small_population")]
    goeaobj = GOEnrichmentStudy(popul_ids, assoc, obo_dag, methods=[method])
    study_ids = [line.rstrip() for line in open(ROOT + "small_study")]
    goea_results = goeaobj.run_study(study_ids, methods=[method])
    return obo_dag, goea_results

if __name__ == '__main__':
    test_zip_nt_lists()
