#!/usr/bin/env python
"""Command-line test for notebooks/relatonships_in_the_go.ipynb."""

from __future__ import print_function

import os
from collections import defaultdict
from goatools.base import get_godag

# pylint: disable=line-too-long
def test_ipyrelatonships_in_the_go():
    """Command-line test for notebooks/relatonships_in_the_go.ipynb."""
    # Loading GO graph with the relationship tags
    repo = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..")
    godag = get_godag(os.path.join(repo, "go-basic.obo"), optional_attrs=['relationship'])

    print("\n## Viewing relationships in the GO graph")
    eg_term = godag['GO:1901990']
    print(set([eg_term]))
    # relationship_rev: 0 items
    # name:regulation of mitotic cell cycle phase transition
    # relationship: 1 items
    #   regulates: 1 items
    #     GO:0044772        level-04        depth-04        mitotic cell cycle phase transition [biological_process]
    # level:6
    # is_obsolete:False
    # namespace:biological_process
    # id:GO:1901990
    # reldepth:7
    # depth:7
    # parents: 2 items
    #   GO:1901987  level-06        depth-06        regulation of cell cycle phase transition [biological_process]
    #   GO:0007346  level-05        depth-05        regulation of mitotic cell cycle [biological_process]
    # children: 6 items
    #   GO:0010389  level-07        depth-08        regulation of G2/M transition of mitotic cell cycle [biological_process]
    #   GO:2000045  level-07        depth-08        regulation of G1/S transition of mitotic cell cycle [biological_process]
    #   GO:0007096  level-07        depth-08        regulation of exit from mitosis [biological_process]
    #   GO:0030071  level-07        depth-10        regulation of mitotic metaphase/anaphase transition [biological_process]
    #   GO:1901991  level-07        depth-08        negative regulation of mitotic cell cycle phase transition [biological_process]
    #   GO:1901992  level-07        depth-08        positive regulation of mitotic cell cycle phase transition [biological_process]
    # _parents: 2 items
    #   GO:1901987
    #   GO:0007346
    # alt_ids: 0 items])


    print("\n## Relationship dictionary")
    print(eg_term.relationship.keys())
    print(eg_term.relationship['regulates'])

    ## Example use case: GO:0007124
    #
    # $ prt_terms GO:0007124
    # [Term]
    # id: GO:0007124
    # name: pseudohyphal growth
    # namespace: biological_process
    # def: "A pattern of cell growth that occurs in conditions of nitrogen limitation
    # and abundant fermentable carbon source. Cells become elongated, switch to a
    # unipolar budding pattern, remain physically attached to each other, and invade
    # the growth substrate." [GOC:krc, PMID:11104818]
    # subset: goslim_candida
    # subset: goslim_yeast
    # is_a: GO:0016049 ! cell growth
    # is_a: GO:0070783 ! growth of unicellular organism as a thread of attached cells
    term_of_interest = godag['GO:0007124']

    # First, find the relationship types which contain "regulates":
    regulates = frozenset([typedef
                           for typedef in godag.typedefs.keys()
                           if 'regulates' in typedef])
    print(regulates)
    assert regulates == frozenset(['regulates', 'negatively_regulates', 'positively_regulates'])

    # Now search through the terms in the tree for those with a relationship in
    # this list and add them to a dictionary dependent on the type of regulation.
    regulating_terms = defaultdict(list)
    for goterm in godag.values():
        if hasattr(goterm, 'relationship'):
            for typedef in regulates.intersection(goterm.relationship.keys()):
                if term_of_interest in goterm.relationship[typedef]:
                    regulating_terms['{:s}d_by'.format(typedef[:-1])].append(goterm)

    # Now regulating_terms contains the GO terms which relate to regulating
    # protein localization to the nucleolus.

    # pseudohyphal growth (GO:0007124) is:
    #
    #   - negatively_regulated_by:
    #     -- GO:2000221 negative regulation of pseudohyphal growth
    #
    #   - regulated_by:
    #     -- GO:2000220 regulation of pseudohyphal growth
    #
    #   - positively_regulated_by:
    #     -- GO:2000222 positive regulation of pseudohyphal growth
    #
    print('{:s} ({:s}) is:'.format(term_of_interest.name, term_of_interest.id))
    for reg_desc, goterms in regulating_terms.items():
        print('\n  - {:s}:'.format(reg_desc))
        for goterm in goterms:
            print('    -- {:s} {:s}'.format(goterm.id, goterm.name))
            for gochild in goterm.children:
                print('    -- {:s} {:s}'.format(gochild.id, gochild.name))


# pseudohyphal growth (GO:0007124) is:
#
#   - negatively_regulated_by:
#     -- GO:2000221 negative regulation of pseudohyphal growth
#     -- GO:0100042 negative regulation of pseudohyphal growth by transcription from RNA polymerase II promoter
#     -- GO:1900462 negative regulation of pseudohyphal growth by negative regulation of transcription from RNA polymerase II promoter
#
#   - regulated_by:
#     -- GO:2000220 regulation of pseudohyphal growth
#
#   - positively_regulated_by:
#     -- GO:2000222 positive regulation of pseudohyphal growth
#     -- GO:0100041 positive regulation of pseudohyphal growth by transcription from RNA polymerase II promoter
#     -- GO:1900461 positive regulation of pseudohyphal growth by positive regulation of transcription from RNA polymerase II promoter


if __name__ == '__main__':
    test_ipyrelatonships_in_the_go()
