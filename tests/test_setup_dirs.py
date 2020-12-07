#!/usr/bin/env python3
"""Test that all GOATOOLS package dirs are in the setup.py file"""

# pylint: disable=wrong-import-position
from os import walk
from os.path import join
from os.path import abspath
import sys
sys.argv = [abspath(__file__), '--help']
from setup import NAME      # goatools
from setup import PACKAGES  # modules in goatools
from tests.utils import REPO


def test_setup_dirs():
    """Test that all GOATOOLS package dirs are in the setup.py file"""
    pkgs_setup = set(m for m in PACKAGES if 'test_' not in m)
    pkgs_dirs = _get_pkgmods()
    assert pkgs_dirs.issubset(pkgs_setup), _errmsg(pkgs_setup, pkgs_dirs)
    print('**NOTE: TEST PASSED')


def _errmsg(pkgs_setup, pkgs_dirs):
    """Print the packages which are not found in setup.py"""
    len_name = len(NAME) + 1
    missing = set(m[len_name:] for m in pkgs_dirs.difference(pkgs_setup))
    return '**FATAL: MISSING PACKAGES in setup.py:\n    NAME + ".{P}",'.format(
        P='",\n    NAME + ".'.join(sorted(missing)))

def _get_pkgmods():
    """Get the GOATOOLS package modules by walking the package dirs"""
    pkgs = set()
    len_repo = len(REPO) + 1
    for root, _, _ in walk(join(REPO, NAME)):
        if root[-11:] != '__pycache__':
            pkg = root[len_repo:].replace('/', '.')
            if 'test_' not in pkg:
                pkgs.add(pkg)
    return pkgs


if __name__ == '__main__':
    test_setup_dirs()
