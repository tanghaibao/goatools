import sys
import os
from goatools.obo_parser import GODag

ROOT = os.path.dirname(os.path.abspath(__file__)) + "/data/"


def print_paths(paths, PRT=sys.stdout):
    for path in paths:
        PRT.write('\n')
        for GO in path:
            PRT.write('{}\n'.format(GO))


def chk_results(actual_paths, expected_paths):
    for actual_path in actual_paths:
        # GOTerm -> list of Strings
        actual = [GO.id for GO in actual_path]
        if actual not in expected_paths:
            raise Exception('ACTUAL {} NOT FOUND IN EXPECTED RESULTS\n'.format(actual))


def test_paths_to_top():
    dag = GODag(ROOT + "mini_obo.obo")
    expected_paths = [['GO:0000001', 'GO:0000002', 'GO:0000005', 'GO:0000010'],
                      ['GO:0000001', 'GO:0000003', 'GO:0000005', 'GO:0000010'],
                      ['GO:0000001', 'GO:0000003', 'GO:0000006', 'GO:0000008', 'GO:0000010']]
    actual_paths = dag.paths_to_top("GO:0000010")
    chk_results(actual_paths, expected_paths)
    print_paths(actual_paths)


def test_paths_to_top_relationships():
    """Test that paths_to_top includes optional relationships (e.g. part_of)."""
    # mini_obo_rel.obo:  root <- B <- D -part_of-> C <- root
    # Without relationships: D -> B -> root  (1 path)
    # With part_of: D -> B -> root  AND  D -> C -> root  (2 paths)
    dag = GODag(ROOT + "mini_obo_rel.obo", optional_attrs={'relationship'})

    # Default behaviour: only is_a paths
    expected_isa_only = [['GO:0000001', 'GO:0000002', 'GO:0000004']]
    actual_isa_only = dag.paths_to_top("GO:0000004")
    chk_results(actual_isa_only, expected_isa_only)
    assert len(actual_isa_only) == 1, actual_isa_only

    # With relationships=True: is_a AND all optional relationships
    expected_with_rels = [
        ['GO:0000001', 'GO:0000002', 'GO:0000004'],
        ['GO:0000001', 'GO:0000003', 'GO:0000004'],
    ]
    actual_with_rels = dag.paths_to_top("GO:0000004", relationships=True)
    chk_results(actual_with_rels, expected_with_rels)
    assert len(actual_with_rels) == 2, actual_with_rels

    # With relationships={'part_of'}: same result as True for this DAG
    actual_part_of = dag.paths_to_top("GO:0000004", relationships={'part_of'})
    chk_results(actual_part_of, expected_with_rels)
    assert len(actual_part_of) == 2, actual_part_of

    print_paths(actual_with_rels)
