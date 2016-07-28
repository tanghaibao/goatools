"""Test get_unique_fields function in nt_utils Python module."""

from goatools.nt_utils import get_unique_fields

def test_get_unique_fields():
    """Test get_unique_fields function in nt_utils Python module."""
    # Test 0. All values are unique
    fld_lists = [
        ['a', 'b', 'c', 'd'],
        ['e', 'f', 'g'],
        ['h', 'i']]
    expected = ['a', 'b', 'c', 'd', 'e', 'f', 'g', 'h', 'i']
    run(fld_lists, expected)
    # Test 1. Field lists have some duplicate fields
    fld_lists = [
        ['a', 'b', 'c', 'd'],
        ['a', 'b', 'g'],
        ['a', 'i']]
    expected = ['a', 'b', 'c', 'd', 'g', 'i']
    run(fld_lists, expected)

def run(fld_lists, expected):
    """Create unique field list. Check for PASS/FAIL."""
    actual = get_unique_fields(fld_lists)
    assert actual == expected, "ACTUAL({A}) != EXPECTED({E})".format(
        A=" ".join(actual), E=" ".join(expected))

if __name__ == '__main__':
    test_get_unique_fields()
