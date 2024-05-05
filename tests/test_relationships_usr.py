#!/usr/bin/env python3
"""Test a user providing both valid and unexpected relationships"""

import os

from goatools.base import get_godag
from goatools.godag.relationship_combos import RelationshipCombos

REPO = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..")


def test_relationships_usr():
    """Test a user providing unexpected relationships"""
    # Set up test
    fin_godag = os.path.join(REPO, "go-basic.obo")
    godag_r1 = get_godag(fin_godag, optional_attrs=["relationship"])
    obj_r1 = RelationshipCombos(godag_r1)
    assert obj_r1.get_set(True) == obj_r1.dag_rels
    assert obj_r1.get_set(False) == set()
    assert obj_r1.get_set(
        set(
            [
                "part_of",
            ]
        )
    ) == {
        "part_of",
    }
    assert obj_r1.get_set(
        {
            "part_of",
        }
    ) == {
        "part_of",
    }
    assert obj_r1.get_set(
        [
            "part_of",
        ]
    ) == {
        "part_of",
    }
    assert obj_r1.get_set("part_of") == {
        "part_of",
    }
    # Run tests: bad user relationships
    assert obj_r1.get_set("BAD_REL") == set()
    assert (
        obj_r1.get_set(
            [
                "BAD_REL",
            ]
        )
        == set()
    )
    assert (
        obj_r1.get_set(
            {
                "BAD_REL",
            }
        )
        == set()
    )
    assert (
        obj_r1.get_set(
            set(
                [
                    "BAD_REL",
                ]
            )
        )
        == set()
    )
    assert obj_r1.get_set(
        set(
            [
                "part_of",
                "BAD_REL",
            ]
        )
    ) == {
        "part_of",
    }

    # Run tests: expected relationships
    godag_r0 = get_godag(fin_godag)
    obj_r0 = RelationshipCombos(godag_r0)
    assert obj_r0.get_set(True) == set()
    assert obj_r0.get_set(False) == set()
    assert (
        obj_r0.get_set(
            set(
                [
                    "part_of",
                ]
            )
        )
        == set()
    )
    assert (
        obj_r0.get_set(
            {
                "part_of",
            }
        )
        == set()
    )
    assert (
        obj_r0.get_set(
            [
                "part_of",
            ]
        )
        == set()
    )
    assert obj_r0.get_set("part_of") == set()
    # Run tests: bad user relationships
    assert obj_r0.get_set("BAD_REL") == set()
    assert (
        obj_r0.get_set(
            [
                "BAD_REL",
            ]
        )
        == set()
    )
    assert (
        obj_r0.get_set(
            {
                "BAD_REL",
            }
        )
        == set()
    )
    assert (
        obj_r0.get_set(
            set(
                [
                    "BAD_REL",
                ]
            )
        )
        == set()
    )
    assert (
        obj_r0.get_set(
            set(
                [
                    "part_of",
                    "BAD_REL",
                ]
            )
        )
        == set()
    )


if __name__ == "__main__":
    test_relationships_usr()
