"""Ensure that alternate GO IDs are in the go-basic.obo DAG go2obj dictionary."""

from goatools.base import get_godag

def test_alt_id():
    """Ensure that alternate GO IDs."""
    obo_dag = get_godag("go-basic.obo", loading_bar=None)
    alt_ids = get_altids(obo_dag)
    obo_goids = obo_dag.keys()
    obo_goids_set = set(obo_goids)
    assert len(alt_ids.intersection(obo_goids_set)) == len(alt_ids)

def get_altids(obo_dag):
    """Get all alternate GO ids for entire go-basic.obo DAG."""
    alt_ids_all = set()
    for _, goobj in obo_dag.items():
        alt_ids_cur = goobj.alt_ids
        if alt_ids_cur:
            alt_ids_all |= set(alt_ids_cur)
    return alt_ids_all

if __name__ == '__main__':
    test_alt_id()
