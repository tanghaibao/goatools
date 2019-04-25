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

import collections as cx

__copyright__ = "Copyright (C) 2016-2019, DV Klopfenstein, H Tang. All rights reserved."
__author__ = "DV Klopfenstein"


# pylint: disable=too-few-public-methods
class AnnotationExtensions(object):
    """A collection of annotation extensions for one gene product."""

    def __init__(self, extensions):
        self.exts = extensions

    def __str__(self):
        hdr = "Ext({N}:{L})".format(
            N=len(self.exts),
            L=",".join(["{}".format(len(ext_lst)) for ext_lst in self.exts]))
        txt = [hdr]
        for ext_lst in self.exts:
            exts_str = ", ".join(str(e) for e in ext_lst)
            txt.append("[{TXT}]".format(TXT=exts_str))
        return " ".join(txt)

    def get_relations_cnt(self):
        """Get the set of all relations."""
        return cx.Counter([e.relation for es in self.exts for e in es])


# Copyright (C) 2016-2019, DV Klopfenstein, H Tang. All rights reserved."
