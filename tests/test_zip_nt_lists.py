"""Test zipping lists whose elements are namedtuples or class objects."""

import os
import collections as cx
from goatools.wr_tbl import zip_nt_lists
from goatools.obo_parser import GODag
from goatools.go_enrichment import GOEnrichmentStudy, get_goea_nts
from goatools.associations import read_associations

def test_zip_nt_lists():
    """Test combining lists whose elements are namedtuples or class objects."""
    ntobj = cx.namedtuple("Nt", "idx")
    goea_results = get_goea_results()
    # Zip a list of namedtuples and another list of namedtuples
    goea_nts = get_goea_nts(goea_results)
    lst2_nts = [ntobj._make([i]) for i in range(len(goea_nts))]
    # Combine lists into a single list whose elements are a namedtuple
    lst_all = zip_nt_lists([lst2_nts, goea_nts])
    assert lst_all[0]._fields == lst2_nts[0]._fields + goea_nts[0]._fields
    # Combine list contains a subset of namedtuple fields
    hdrs = ['idx', 'NS', 'level', 'depth', 'GO',
            'study_count', 'study_n', 'pop_count', 'pop_n', 'p_fdr_bh', 'name']
    lst_sub = zip_nt_lists([lst2_nts, goea_nts], hdrs)
    assert list(lst_sub[0]._fields) == hdrs, "{F} {H}".format(F=lst_sub[0]._fields, H=hdrs)


def get_goea_results(method="fdr_bh"):
    """Get GOEA results."""
    root_dir = os.path.dirname(os.path.abspath(__file__)) + "/data/"
    obo_dag = GODag(root_dir + "goslim_generic.obo")
    assoc = read_associations(root_dir + "slim_association", no_top=True)
    popul_ids = [line.rstrip() for line in open(root_dir + "small_population")]
    goeaobj = GOEnrichmentStudy(popul_ids, assoc, obo_dag, methods=[method])
    study_ids = [line.rstrip() for line in open(root_dir + "small_study")]
    goea_results = goeaobj.run_study(study_ids, methods=[method])
    return goea_results

if __name__ == '__main__':
    test_zip_nt_lists()
