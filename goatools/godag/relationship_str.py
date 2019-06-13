"""Create strings representing relationships on GO Terms.

                   +------- has 'part_of' relationship(s)
                   |    +-- pointed to by a GO ID with a 'part_of' relationship
                   |    |
                   V    V
GO:0008150 L00 D00 .... .rdu biological_process
GO:0050896 L01 D01 .... .rdu response to stimulus
GO:0042221 L02 D02 .... p... response to chemical
GO:0032501 L01 D01 .... .rdu multicellular organismal process
GO:0003008 L02 D02 .... .r.. system process
GO:0051606 L02 D02 .... .... detection of stimulus
GO:0050877 L03 D03 .... .rdu nervous system process
GO:0009593 L03 D03 P... .... detection of chemical stimulus
GO:0007600 L04 D04 .... pr.. sensory perception
GO:0050906 L03 D03 P... .... detection of stimulus involved in sensory perception
GO:0050890 L04 D04 .... .... cognition
GO:0050907 L04 D04 P... .... detection of chemical stimulus involved in sensory perception
GO:0007606 L05 D05 .... p... sensory perception of chemical stimulus
GO:0050893 L05 D05 P... .... sensory processing
GO:0050911 L05 D05 P... .... detection of chemical stimulus involved in sensory perception of smell
GO:0007608 L06 D06 .... p... sensory perception of smell

"""

__copyright__ = "Copyright (C) 2010-2019, DV Klopfenstein, H Tang, All rights reserved."
__author__ = "DV Klopfenstein"

from collections import OrderedDict
from goatools.godag.consts import RELATIONSHIP_LIST
from goatools.godag.consts import RELATIONSHIP_SET

# pylint: disable=too-few-public-methods,bad-whitespace
class RelationshipStr(object):
    """Create strings representing relationships on GO Terms."""

    # go-basic.obo: fmt(1.2) rel(2019-02-20) 47,177 GO Terms; optional_attrs(relationship)
    # relationship:
    #     6,882 part_of
    #     3,230 regulates
    #     2,804 negatively_regulates
    #     2,785 positively_regulates

    rel2chr = OrderedDict([
        ('part_of',   'P'),
        ('regulates', 'R'),
        ('negatively_regulates', 'D'),
        ('positively_regulates', 'U')])

    rev2chr = OrderedDict([
        ('part_of',   'p'),
        ('regulates', 'r'),
        ('negatively_regulates', 'd'),
        ('positively_regulates', 'u')])

    def __init__(self, relationships=None):
        assert set(self.rel2chr.keys()) == RELATIONSHIP_SET
        # Ordered relationships
        _rels = relationships if relationships else set()
        self.rels = [r for r in RELATIONSHIP_LIST if r in _rels]

    def str_relationships(self, goobj):
        """Get a string representing the presence of absence of relationships. Ex: P..."""
        rel_cur = goobj.relationship
        return "".join([self.rel2chr.get(r, '?') if r in rel_cur else '.' for r in self.rels])

    def str_rel_short(self, goobj):
        """Get a string representing the presence of absence of relationships. Ex: P"""
        if not goobj.relationship:
            return ''
        rel_cur = goobj.relationship
        return "".join([self.rel2chr.get(r, '?') for r in self.rels if r in rel_cur])

    def str_relationships_rev(self, goobj):
        """Get a string representing the presence of absence of relationships. Ex: pr.."""
        rel_cur = goobj.relationship_rev
        return "".join([self.rev2chr[r] if r in rel_cur else '.' for r in self.rels])

    def prt_keys(self, prt, pre):
        """Print the alias for a relationship and its alias."""
        prt.write('{PRE}Relationship to parent: {ABC}\n'.format(
            PRE=pre, ABC=''.join(self.rel2chr.values())))
        for rel, alias in self.rel2chr.items():
            prt.write('{PRE}    {A} {DESC}\n'.format(PRE=pre, A=alias, DESC=rel))
        prt.write('\n{PRE}Relationship to child: {ABC}\n'.format(
            PRE=pre, ABC=''.join(self.rev2chr.values())))
        for rel, alias in self.rev2chr.items():
            prt.write('{PRE}    {A} {DESC}\n'.format(PRE=pre, A=alias, DESC=rel))


# Copyright (C) 2010-2019, DV Klopfenstein, H Tang, All rights reserved.
