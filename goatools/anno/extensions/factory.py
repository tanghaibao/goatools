"""Annotation Extension for relational expressions.

   https://link.springer.com/protocol/10.1007/978-1-4939-3743-1_17

   Correlated function associations between genes and GOs containing contextual information.

   With contextual information identify gene products that perform a role:
    - only under certain conditions
    - in the presence of specific factors

   Gene products can have different roles:
     - in different cells
     - in different tissues

"""

import sys
from goatools.anno.extensions.extensions import AnnotationExtensions
from goatools.anno.extensions.extension import AnnotationExtension

__copyright__ = "Copyright (C) 2016-2019, DV Klopfenstein, H Tang. All rights reserved."
__author__ = "DV Klopfenstein"


def get_extensions(extstr):
    """Return zero or greater Annotation Extensions, given a line of text."""
    # Extension examples:
    #   has_direct_input(UniProtKB:P37840),occurs_in(GO:0005576)
    #   part_of(UBERON:0006618),part_of(UBERON:0002302)
    #   occurs_in(CL:0000988)|occurs_in(CL:0001021)
    if not extstr:
        return None
    exts = []
    for ext_lst in extstr.split('|'):
        grp = []
        for ext in ext_lst.split(','):
            idx = ext.find('(')
            if idx != -1 and ext[-1] == ')':
                grp.append(AnnotationExtension(ext[:idx], ext[idx+1:-1]))
            else:
                # Ignore improperly formatted Extensions
                sys.stdout.write('BAD Extension({E})\n'.format(E=ext))
        exts.append(grp)
    return AnnotationExtensions(exts)


# Copyright (C) 2016-2019, DV Klopfenstein, H Tang. All rights reserved."
