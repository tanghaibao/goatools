"""Annotation Extension for relational expressions.

   https://link.springer.com/protocol/10.1007/978-1-4939-3743-1_17

"""

__copyright__ = "Copyright (C) 2016-2019, DV Klopfenstein, H Tang. All rights reserved."
__author__ = "DV Klopfenstein"

# pylint: disable=too-few-public-methods
class AnnotationExtension(object):
    """Annotation Extension for relational expressions."""

    def __init__(self, relation, entity):
        # Relationship between GO term and the entity
        self.relation = relation
        # An identifier for a database object or ontology term
        self.entity = entity

    def __str__(self):
        return "{R}({E})".format(R=self.relation, E=self.entity)


# Copyright (C) 2016-2019, DV Klopfenstein, H Tang. All rights reserved."
