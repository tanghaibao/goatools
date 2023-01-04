#!/usr/bin/env python
"""Test Yang's RWC measure added to other semantic similariy measures"""


from os.path import join
from os.path import dirname
from os.path import abspath
from goatools.obo_parser import GODag
from goatools.semantic import deepest_common_ancestor
from goatools.semantic import semantic_distance


REPO = join(dirname(abspath(__file__)), "..")


def test_i256_mca_dist():
    """Test faster version of sematic similarity"""
    godag = GODag(join(REPO, 'tests/data/yangRWC/fig1b.obo'))

    go_ids = [
        'GO:0000006',
        'GO:0000007',
    ]

    mca = deepest_common_ancestor(go_ids, godag)
    dist = semantic_distance(*go_ids, godag)

    print(f'{mca} is the Most Recent Common Ancestor of {go_ids}')
    print(f'{dist} minimum number of connecting branches aka semantic distance between {go_ids}')

    assert mca == 'GO:0000014'
    assert dist == 4


if __name__ == '__main__':
    test_i256_mca_dist()
