"""Plot a GODagSmall."""

__copyright__ = "Copyright (C) 2016, DV Klopfenstein, H Tang, All rights reserved."
__author__ = "DV Klopfenstein"

import sys 
from collections import OrderedDict
from goatools.wr_tbl import get_fmtflds

class GODagSmallPlot(object):
    """Plot a GODagSmall."""

    # http://www.graphviz.org/doc/info/colors.html
    rel2col = {
        'is_a':      'black',
        'part_of':   'blue',
        'regulates': 'gold',
        'positively_regulates': 'green',
        'negatively_regulates': 'red',
        'occurs_in':            'aquamarine4',
        'capable_of':           'dodgerblue',
        'capable_of_part_of':   'darkorange',
    }

    alpha2col = OrderedDict([
        (0.005, 'Yellow4'),
        (0.01,  'Yellow3'),
        (0.05,  'Yellow2'),
    ])

    fmtstr = "{GO} L{level:>02} D{depth:>02}\n{name}"

    def __init__(self, godagsmall, **kws):
        # init fmtstr_usr
        self.fmtstr_usr = None
        self.fmtfld_usr = None
        self._init_fmtstr_usr(**kws)
        #
        self.log = kws['log'] if 'log' in kws else sys.stdout
        self.go2nt = kws['go2nt'] if 'go2nt' in kws else None
        self.alpha_str = kws['alpha_str'] if 'alpha_str' in kws else None
        self.godag = godagsmall
        self.dpi = 150
        self.pydot = None

    def _init_fmtstr_usr(self, **kws):
        if 'fmtstr' in kws:
            self.fmtstr_usr = kws['fmtstr']
            self.fmtfld_usr = get_fmtflds(self.fmtstr_usr)
        return None

    # ----------------------------------------------------------------------------------
    # pydot
    def plt_pydot(self, fout_png):
        dag = self.get_pydot_graph()
        dag.write_png(fout_png)
        self.log.write("  WROTE: {F}\n".format(F=fout_png))

    def get_pydot_graph(self):
        rel = "is_a"
        pydot = self.get_pydot()
        # Initialize empty dag
        dag = pydot.Dot(graph_type='digraph', dpi="{}".format(self.dpi)) # Directed Graph
        # Initialize nodes
        go2node = self.get_go2pydotnode()
        # Add nodes to graph
        for node in go2node.values():
            dag.add_node(node)
        # Add edges to graph
        for src, tgt in self.godag.get_edges():
            dag.add_edge(pydot.Edge(go2node[tgt], go2node[src],
                shape = "normal",
                color = self.rel2col[rel],
                dir = "back")) # invert arrow direction for obo dag convention
        return dag

    def get_go2pydotnode(self):
        """Get pydot Nodes.""" 
        go2node = {}
        for goid, goobj in self.godag.go2obj.items():
            txt = self.get_node_text(goid, goobj)
            fillcolor, colorscheme = self._get_fillcolor(goid)
            node = self.pydot.Node(
                txt,
                shape = "box",
                style = "rounded, filled",
                fillcolor = fillcolor,
                #colorscheme = colorscheme,
                color = "mediumseagreen")
            go2node[goid] = node
        return go2node

    def _get_fillcolor(self, goid):
        # fillcolor default
        fillcolor = "white"
        colorscheme = None
        # fillcolor based on pvalue
        if self.alpha_str is not None and self.go2nt is not None and goid in self.go2nt:
            nt = self.go2nt[goid]
            pval = getattr(nt, self.alpha_str, None)
            if pval is not None:
                for alpha, color in self.alpha2col.items():
                    if pval < alpha:
                        return color, colorscheme
        return fillcolor, colorscheme

    def get_node_text(self, goid, goobj):
        #if self.fmtfld_usr is not None and self.go2nt is not None and goid in self.go2nt:
        txt = self.fmtstr.format(
            GO = goobj.id.replace("GO:", ""),
            level = goobj.level,
            depth = goobj.depth,
            name = goobj.name)
        return txt.replace(",", "\n")

    def get_pydot(self):
        """Return pydot package. Load pydot, if necessary."""
        if self.pydot:
            return self.pydot
        self.pydot = __import__("pydot")
        return self.pydot

# Copyright (C) 2016, DV Klopfenstein, H Tang, All rights reserved.
