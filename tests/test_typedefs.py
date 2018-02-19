#!/usr/bin/env python
"""Ensure that alternate GO IDs are in the go-basic.obo DAG go2obj dictionary."""

from goatools.base import get_godag


def test_typedef():
    """Ensure that alternate GO IDs."""
    obo_dag = get_godag("go-basic.obo", loading_bar=None)
    print(obo_dag.typedefs['negatively_regulates'])


if __name__ == '__main__':
    test_typedef()
