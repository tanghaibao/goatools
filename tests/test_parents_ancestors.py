#!/usr/bin/env python
"""Test initializing TermCounts with annotations made to alternate GO ID"""


import os
from goatools.obo_parser import GODag
# from goatools.anno.idtogos_reader import IdToGosReader
from goatools.gosubdag.gosubdag import GoSubDag
from goatools.godag.go_tasks import get_go2children
from goatools.godag.go_tasks import get_go2parents


REPO = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..")

def test_parents_ancestors():
    """Test getting parents and ancestors"""
    # Load a small GO DAG to demonstrate getting parents and ancestors
    file_dag = os.path.join(REPO, 'tests/data/i126/viral_gene_silence.obo')
    # Load all relationships using optional attribute
    godag = GODag(file_dag)

    optional_relationships = set()  # Don't trace any optional relationships
    go2parents_isa = get_go2parents(godag, optional_relationships)
    go2children_isa = get_go2children(godag, optional_relationships)
    # TODO: Add more tests for only is_a


    godag = GODag(file_dag, optional_attrs={'relationship'})
    goids = set(o.item_id for o in godag.values())

    # Get parents through "is_a" only
    optional_relationships = set()  # Don't trace any optional relationships
    go2parents_isa = get_go2parents(godag, optional_relationships)
    go2children_isa = get_go2children(godag, optional_relationships)

    # Get parents through "is_a" and all the "regulates" realtionships
    optional_relationships = {'regulates', 'negatively_regulates', 'positively_regulates'}
    go2parents_reg = get_go2parents(godag, optional_relationships)
    go2children_reg = get_go2children(godag, optional_relationships)

    # Print parents throush "is_a" relationship
    goid = 'GO:0019222'  # regulation of metabolic process
    assert go2parents_isa[goid] == {'GO:0050789'}
    assert go2parents_reg[goid] == {'GO:0050789', 'GO:0008152'}

    exp = {'GO:0009892', 'GO:0060255'}
    assert go2children_isa[goid] == exp
    assert go2children_reg[goid] == exp
    assert go2children_isa['GO:0008152'] == {'GO:0071704'}
    assert go2children_reg['GO:0008152'] == {'GO:0071704', 'GO:0019222', 'GO:0009892'}

    # Load GO DAG into a GoSubDag object, to use user-selected relationships
    gosubdag_r0 = GoSubDag(goids, godag)
    assert gosubdag_r0.rcntobj.go2ancestors[goid] == \
        {'GO:0050789', 'GO:0065007', 'GO:0008150'}

    # Load GO DAG into a GoSubDag object, to use user-selected relationships
    gosubdag_r1 = GoSubDag(goids, godag, relationships=optional_relationships)
    assert gosubdag_r1.rcntobj.go2ancestors[goid] == \
        {'GO:0050789', 'GO:0008152', 'GO:0065007', 'GO:0008150'}, \
        gosubdag_r1.rcntobj.go2ancestors[goid]

    exp = {'GO:0071704', 'GO:0010467', 'GO:0043170'}
    assert gosubdag_r0.rcntobj.go2descendants['GO:0008152'] == exp

    assert gosubdag_r0.rcntobj.go2descendants['GO:0043170'] == {'GO:0010467'}
    exp = {'GO:0010467', 'GO:0010468', 'GO:0010605', 'GO:0010608', 'GO:0010629',
           'GO:0016441', 'GO:0016458', 'GO:0040029', 'GO:0060147', 'GO:0060148',
           'GO:0060150', 'GO:0060255', 'GO:0060968'}
    assert gosubdag_r1.rcntobj.go2descendants['GO:0043170'] == exp

    gosubdag_r1n = GoSubDag(goids, godag, relationships={'negatively_regulates'})
    exp = {'GO:0010629', 'GO:0016441', 'GO:0016458'}
    assert gosubdag_r1n.rcntobj.go2descendants['GO:0010467'] == exp

    ## # Load annotation data
    ## file_anno = os.path.join(REPO, './notebooks/data/sorghumannotation.anno')
    ## annoobj = IdToGosReader(file_anno, godag=godag)
    ## # Get annotations: geneid-to-set_of_GO_IDs
    ## id2gos = annoobj.get_id2gos('BP')


if __name__ == '__main__':
    test_parents_ancestors()
