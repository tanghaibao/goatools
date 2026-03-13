#!/usr/bin/env python3
"""Test running go_plot subdag with different output orientation."""

from __future__ import print_function

import os
import subprocess
import tempfile


def cmds_plot_annos():
    """Test running go_plot subdag with different output orientation."""
    with tempfile.TemporaryDirectory() as tmpdir:
        for idx, (cmd, check) in enumerate(_get_cmds(tmpdir)):
            print('------------------- TEST {I} ------------------------------------'.format(I=idx))
            print(f'CMD: {cmd}')
            r = subprocess.run(cmd, shell=True, cwd='..')
            assert r.returncode == 0

            with open(os.path.join(tmpdir, 'bbbb.dot')) as f:
                ff = f.read()
                if check not in ff:
                    raise AssertionError(f'Expected phrase "{check}" were not found in dot output.')

            print("TEST PASSED")


def _get_cmds(tmpdir):
    """Commands for generation of different styles, output md5sum"""
    # pylint: disable=line-too-long
    fout_dot = os.path.join(tmpdir, 'bbbb.dot')
    return [
        ('python3 -m goatools go_plot GO:0000010 -o {DOT} --obo=tests/data/mini_obo.obo'.format(DOT=fout_dot), 'rankdir=TB'),
        ('python3 -m goatools go_plot GO:0000010 -o {DOT} --obo=tests/data/mini_obo.obo --rankdir=TB'.format(DOT=fout_dot), 'rankdir=TB'),
        ('python3 -m goatools go_plot GO:0000010 -o {DOT} --obo=tests/data/mini_obo.obo --rankdir=LR'.format(DOT=fout_dot), 'rankdir=LR'),
        ('python3 -m goatools go_plot GO:0000010 -o {DOT} --obo=tests/data/mini_obo.obo --rankdir=BT'.format(DOT=fout_dot), 'rankdir=BT'),
        ('python3 -m goatools go_plot GO:0000010 -o {DOT} --obo=tests/data/mini_obo.obo --rankdir=RL'.format(DOT=fout_dot), 'rankdir=RL'),
    ]


if __name__ == '__main__':
    cmds_plot_annos()
