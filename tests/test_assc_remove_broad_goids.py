#!/usr/bin/env python
"""Test choosing a set of broad GO IDs to be removed"""

from goatools.anno.update_association import get_goids_to_remove
from goatools.anno.update_association import remove_assc_goids
from goatools.anno.broad_gos import NS2GOS_SHORT
from goatools.anno.broad_gos import NS2GOS
from goatools.godag.consts import TOP_TERMS


def test_remove_assc_goids():
    """Test choosing a set of broad GO IDs to be removed"""
    # The default is to remove a small number of broad GO IDs
    assert get_goids_to_remove() == set.union(*NS2GOS_SHORT.values())
    assert get_goids_to_remove(None) == set.union(*NS2GOS_SHORT.values())
    # Remove a larger number of borad GO IDs (~80 out of 45k)
    assert get_goids_to_remove(True) == set.union(*NS2GOS.values())
    # Do not remove any GO IDs
    assert get_goids_to_remove(False) == set()
    # Remove user-specified GO IDs
    assert get_goids_to_remove(TOP_TERMS) == TOP_TERMS
    assert get_goids_to_remove(list(TOP_TERMS)) == TOP_TERMS

    assoc = {1:{'a', 'b', 'c'}, 2:{'b', 'c', 'd'}, 3:{'a', 'b'}}

    # pylint: disable=line-too-long
    ret = remove_assc_goids(assoc, {'a'})
    assert ret == {'goids_removed': set(['a']),
                   'genes_removed': set([]),
                   'assoc_reduced': {1: set(['c', 'b']), 2: set(['c', 'b', 'd']), 3: set(['b'])}}, ret

    ret = remove_assc_goids(assoc, {'a', 'b'})
    assert ret == {'goids_removed': set(['a', 'b']),
                   'genes_removed': set([3]),
                   'assoc_reduced': {1: set(['c']), 2: set(['c', 'd'])}}, ret

    ret = remove_assc_goids(assoc, {'a', 'b', 'x'})
    assert ret == {'goids_removed': set(['a', 'b']),
                   'genes_removed': set([3]),
                   'assoc_reduced': {1: set(['c']), 2: set(['c', 'd'])}}


if __name__ == '__main__':
    test_remove_assc_goids()
