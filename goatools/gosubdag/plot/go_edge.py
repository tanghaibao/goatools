"""Create a pydot Edge for a GO Term."""

__copyright__ = "Copyright (C) 2019-2020, DV Klopfenstein, H Tang, All rights reserved."
__author__ = "DV Klopfenstein"

import itertools
import pydot
from goatools.gosubdag.utils import extract_kwargs


#pylint: disable=too-few-public-methods
class GoEdgeOpts:
    """Processes edges in GO plot args."""

    exp_keys = set(['edge2txt',])
    exp_elems = None  # set()

    def __init__(self, gosubdag, **kws):
        self.gosubdag = gosubdag
        # kws = {'set':set(...), 'dict':{...}}
        self.kws = extract_kwargs(kws, self.exp_keys, self.exp_elems)

    def get_kws(self):
        """Only load keywords if they are specified by the user."""
        ret = self.kws['dict'].copy()
        if 'dict' in self.kws and 'edge2txt' in self.kws['dict']:
            self._init_edge2txt_altgos(self.kws['dict']['edge2txt'])
        return ret

    def _init_edge2txt_altgos(self, edge2txt):
        """If user provided GO.alt_id, add the corressponding main GO ID, if needed"""
        _go2obj = self.gosubdag.go2obj
        for edge, txt in edge2txt.items():
            go_lo_usr, go_hi_usr = edge
            go_lo_hi_sets = [
                set([go_lo_usr, _go2obj[go_lo_usr].item_id]),
                set([go_hi_usr, _go2obj[go_hi_usr].item_id])]
            for go_lo_hi_cur in itertools.product(*go_lo_hi_sets):
                if go_lo_hi_cur not in edge2txt:
                    ## print('ADDING ALT EDGE:', go_lo_hi_cur)
                    edge2txt[edge] = txt


class GoEdge:
    """Creates pydot Edge containing a GO term."""

    def __init__(self, gosubdag, optobj):
        self.gosubdag = gosubdag     # GoSubDag
        self.kws = optobj.get_kws()  # GoEdgeOpts -> text  options
    ##     self.present = optobj.get_present()
    ##     self.objcolor = objcolor     # Go2Color   -> color options
    ##     self.go2color = objcolor.go2color

    def add_edges(self, edges_list, go2node, dag, **kws):
        """Add Edges from one Node to another Node to Dag"""
        # style: solid dashed dotted bold invis tapered
        # arrowType: http://www.graphviz.org/doc/info/attrs.html#k:arrowType
        edge2txt = self.kws.get('edge2txt')
        for src, tgt in edges_list:
            assert src in go2node, "MISSING Edge source({S}); target({T})".format(S=src, T=tgt)
            assert tgt in go2node, "MISSING Edge target({T}); source({S})".format(S=src, T=tgt)
            if edge2txt and (src, tgt) in edge2txt:
                kws = dict(kws)
                kws['label'] = edge2txt[(src, tgt)]
            dag_edge = pydot.Edge(
                go2node[tgt], go2node[src],
                shape="normal",
                # # style="normal",
                # color=color,
                dir="back", # invert arrow direction for obo dag convention
                **kws)
            # sequence parent_graph points attributes type parent_edge_list
            # GoSubDagPlot._prt_edge(dag_edge, 'parent_edge_list')
            dag.add_edge(dag_edge)


# Copyright (C) 2019-2020, DV Klopfenstein, H Tang, All rights reserved.
