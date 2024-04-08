#!/usr/bin/env python3
"""Test running go_plot subdag with different output orientation."""

from __future__ import print_function

import os
import subprocess


def cmds_plot_annos():
    """Test running go_plot subdag with different output orientation."""
    for idx, (cmd, check) in enumerate(_get_cmds()):
        print('------------------- TEST {I} ------------------------------------'.format(I=idx))
        print(f'CMD: {cmd}')
        r = subprocess.run(cmd, shell=True, cwd='..')
        assert r.returncode == 0

        with open('bbbb.dot') as f:
            ff = f.read()
            if check not in ff:
                raise AssertionError(f'Expected phrase "{check}" were not found in dot output.')
        os.remove('bbbb.dot')

        print("TEST PASSED")


def _get_cmds():
    """Commands for generation of different styles, output md5sum"""
    # pylint: disable=line-too-long
    return [
        ('scripts/go_plot.py GO:0000010 -o tests/bbbb.dot --obo=tests/data/mini_obo.obo', 'rankdir=TB'),
        ('scripts/go_plot.py GO:0000010 -o tests/bbbb.dot --obo=tests/data/mini_obo.obo --rankdir=TB', 'rankdir=TB'),
        ('scripts/go_plot.py GO:0000010 -o tests/bbbb.dot --obo=tests/data/mini_obo.obo --rankdir=LR', 'rankdir=LR'),
        ('scripts/go_plot.py GO:0000010 -o tests/bbbb.dot --obo=tests/data/mini_obo.obo --rankdir=BT', 'rankdir=BT'),
        ('scripts/go_plot.py GO:0000010 -o tests/bbbb.dot --obo=tests/data/mini_obo.obo --rankdir=RL', 'rankdir=RL'),
    ]


if __name__ == '__main__':
    cmds_plot_annos()
