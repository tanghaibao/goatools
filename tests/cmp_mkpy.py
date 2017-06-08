#!/usr/bin/env python
"""Determine which tests are in the directory which should be added to the makefile."""

import os
import sys
import re
import glob

def main():
    """Compare the tests in this directory to the tests in the makefile."""
    tests_cwd = set(glob.glob('*.py'))
    tests_mk = _get_makefile_tests()
    _prt_tests_to_add_to_mk(tests_cwd, tests_mk)

def _prt_tests_to_add_to_mk(tests_cwd, tests_mk):
    """Print tests to add to the makefile."""
    tests = tests_cwd.symmetric_difference(tests_mk)
    exclude = set([__file__[2:], "test_goea_errors.py"])
    if tests:
        sys.stdout.write("ADD THESE TESTS TO THE MAKEFILE:\n")
        for test in tests:
            if test not in exclude:
                sys.stdout.write("\t$(PYTHON) {TEST}\n".format(TEST=test))

def _get_makefile_tests():
    """Get the list of tests in the makefile."""
    fin_makefile = "makefile"
    cwd = os.path.dirname(os.path.abspath(__file__))
    with open(os.path.join(cwd, fin_makefile)) as ifstrm:
        tests = set()
        for line in ifstrm:
            if line[:11] == "\t$(PYTHON) ":
                mtch = re.match(r'(\S+.py)', line[11:])
                if mtch:
                    tests.add(mtch.group(1))
        sys.stdout.write("READ: {MK}\n".format(MK=fin_makefile))
    return tests

if __name__ == '__main__':
    main()
