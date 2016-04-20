"""Plot a GODagSmall."""

__copyright__ = "Copyright (C) 2016, DV Klopfenstein, H Tang, All rights reserved."
__author__ = "DV Klopfenstein"

import sys 
import os
import collections as cx
from collections import OrderedDict
from goatools.wr_tbl import get_fmtflds
from goatools.godag_obosm import OboToGoDagSmall

def plot_gos(fout_png, goids, obo_dag, *args, **kws):
    """Given GO ids and the obo_dag, create a plot of paths from GO ids."""
    engine = kws['engine'] if 'engine' in kws else 'pydot'
    godagsmall = OboToGoDagSmall(goids=goids, obodag=obo_dag).godag
    godagplot = GODagSmallPlot(godagsmall, *args, **kws)
    godagplot.plt(fout_png, engine)

def plot_goid2goobj(fout_png, goid2goobj, *args, **kws):
    """Given a dict containing GO id and its goobj, create a plot of paths from GO ids."""
    engine = kws['engine'] if 'engine' in kws else 'pydot'
    godagsmall = OboToGoDagSmall(goid2goobj=goid2goobj).godag
    godagplot = GODagSmallPlot(godagsmall, *args, **kws)
    godagplot.plt(fout_png, engine)

def plot_results(fout_png, goea_results, *args, **kws):
    """Given a list of GOEA results, plot result GOs up to top."""
    if "{NS}" not in fout_png:
        plt_goea_results(fout_png, goea_results, *args, **kws)
    else:
        # Plot separately by NS: BP, MF, CC
        ns2goea_results = cx.defaultdict(list)
        for rec in goea_results:
            ns2goea_results[rec.NS].append(rec)
        for ns, ns_res in ns2goea_results.items():
            png = fout_png.format(NS=ns)
            plt_goea_results(png, ns_res, *args, **kws)

def plt_goea_results(fout_png, goea_results, *args, **kws):
    """Plot a single page."""
    engine = kws['engine'] if 'engine' in kws else 'pydot'
    godagsmall = OboToGoDagSmall(goea_results=goea_results).godag
    godagplot = GODagSmallPlot(godagsmall, *args, goea_results=goea_results, **kws)
    godagplot.plt(fout_png, engine)

class GODagPltVars(object):

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
        # GOEA GO terms that are significant
        (0.005, 'mistyrose'),
        (0.01,  'moccasin'),
        (0.05,  'lemonchiffon1'),
        # GOEA GO terms that are not significant
        (1.00,  'grey95'),
    ])

    key2col = {
        'level_01': 'lightcyan',
        'go_sources': 'palegreen',
    }

    fmthdr = "{GO} L{level:>02} D{depth:>02}"
    fmtres = "{study_count} genes"
    # study items per line on GO Terms:
    items_p_line = 5 


class GODagSmallPlot(object):
    """Plot a graph contained in an object of type GODagSmall ."""

    def __init__(self, godagsmall, *args, **kws):
        self.args = args
        self.log = kws['log'] if 'log' in kws else sys.stdout
        self.go2res = self._init_go2res(**kws)
        self.id2symbol = kws['id2symbol'] if 'id2symbol' in kws else {}
        self.study_items = kws['study_items'] if 'study_items' in kws else None
        self.study_items_max = self._init_study_items_max()
        self.alpha_str = kws['alpha_str'] if 'alpha_str' in kws else None
        self.pltvars = kws['GODagPltVars'] if 'GODagPltVars' in kws else GODagPltVars()
        if 'items_p_line' in kws:
            self.pltvars.items_p_line = kws['items_p_line']
        self.dpi = kws['dpi'] if 'dpi' in kws else 150
        self.godag = godagsmall
        self.goid2color = self._init_goid2color()
        self.pydot = None

    def _init_study_items_max(self):
        """User can limit the number of genes printed in a GO term."""
        if self.study_items is None:
            return None
        if self.study_items is True:
            return None
        if isinstance(self.study_items, int):
            return self.study_items
        return None

    def _init_go2res(self, **kws):
        if 'goea_results' in kws:
            return {res.GO:res for res in kws['goea_results']}

    def _init_goid2color(self):
        """Set colors of GO terms."""
        goid2color = {}
        # 1. colors based on p-value override colors based on source GO
        alpha2col = self.pltvars.alpha2col
        if self.go2res is not None:
            for res in self.go2res.values():
                pval = res.get_pvalue()
                for alpha, color in alpha2col.items():
                    if pval <= alpha and res.study_count != 0:
                        goid = getattr(res, 'GO')
                        if goid not in goid2color:
                            goid2color[goid] = color
        # 2. GO source color
        color = self.pltvars.key2col['go_sources']
        for goid in self.godag.go_sources:
            if goid not in goid2color:
                goid2color[goid] = color
        # 3. Level-01 GO color
        color = self.pltvars.key2col['level_01']
        for goid, goobj in self.godag.go2obj.items():
            if goobj.level == 1: 
                if goid not in goid2color:
                    goid2color[goid] = color
        return goid2color

    def plt(self, fout_img, engine="pydot"):
        """Plot using pydot, graphviz, or GML."""
        if engine == "pydot":
            self._plt_pydot(fout_img)
        elif engine == "pygraphviz":
            raise Exception("TO BE IMPLEMENTED SOON: ENGINE pygraphvis")
        else:
            raise Exception("UNKNOWN ENGINE({E})".format(E=engine))

    # ----------------------------------------------------------------------------------
    # pydot
    def _plt_pydot(self, fout_img):
        """Plot using the pydot graphics engine."""
        dag = self._get_pydot_graph()
        img_fmt = os.path.splitext(fout_img)[1][1:]
        dag.write(fout_img, format=img_fmt)
        self.log.write("  WROTE: {F}\n".format(F=fout_img))

    def _get_pydot_graph(self):
        """Given a DAG, return a pydot digraph object."""
        rel = "is_a"
        pydot = self._get_pydot()
        # Initialize empty dag
        dag = pydot.Dot(graph_type='digraph', dpi="{}".format(self.dpi)) # Directed Graph
        # Initialize nodes
        go2node = self._get_go2pydotnode()
        # Add nodes to graph
        for node in go2node.values():
            dag.add_node(node)
        # Add edges to graph
        rel2col = self.pltvars.rel2col
        for src, tgt in self.godag.get_edges():
            dag.add_edge(pydot.Edge(go2node[tgt], go2node[src],
                shape = "normal",
                color = rel2col[rel],
                dir = "back")) # invert arrow direction for obo dag convention
        return dag

    def _get_go2pydotnode(self):
        """Create pydot Nodes.""" 
        go2node = {}
        for goid, goobj in self.godag.go2obj.items():
            txt = self._get_node_text(goid, goobj)
            fillcolor = self.goid2color.get(goid, "white")
            node = self.pydot.Node(
                txt,
                shape = "box",
                style = "rounded, filled",
                fillcolor = fillcolor,
                color = "mediumseagreen")
            go2node[goid] = node
        return go2node

    def _get_pydot(self):
        """Return pydot package. Load pydot, if necessary."""
        if self.pydot:
            return self.pydot
        self.pydot = __import__("pydot")
        return self.pydot

    # ----------------------------------------------------------------------------------
    # Methods for text printed inside GO terms
    def _get_node_text(self, goid, goobj):
        """Return a string to be printed in a GO term box."""
        txt = []
        # Header line: "GO:0036464 L04 D06"
        txt.append(self.pltvars.fmthdr.format(
            GO = goobj.id.replace("GO:", "GO"),
            level = goobj.level,
            depth = goobj.depth))
        # GO name line: "cytoplamic ribonucleoprotein"
        name = goobj.name.replace(",", "\n")
        txt.append(name)
        # study info line: "24 genes"
        study_txt = self._get_study_txt(goid)
        if study_txt is not None:
            txt.append(study_txt)
        # return text string
        return "\n".join(txt)

    def _get_study_txt(self, goid):
        """Get GO text from GOEA study."""
        if self.go2res is not None:
            res = self.go2res.get(goid, None)
            if res is not None:
                if self.study_items is not None:
                    return self._get_item_str(res)
                else:
                    return self.pltvars.fmtres.format(
                        study_count = res.study_count)

    def _get_item_str(self, res):
        """Return genes in any of these formats:
              1. 19264, 17319, 12520, 12043, 74131, 22163, 12575
              2. Ptprc, Mif, Cd81, Bcl2, Sash3, Tnfrsf4, Cdkn1a
              3. 7: Ptprc, Mif, Cd81, Bcl2, Sash3...
        """
        N = self.pltvars.items_p_line
        prt_items = sorted([self.__get_genestr(itemid) for itemid in res.study_items])
        prt_multiline = [prt_items[i:i+N] for i in range(0, len(prt_items), N)]
        if self.study_items_max is None:
            return "\n".join([", ".join(str(e) for e in sublist) for sublist in prt_multiline])
        else:
            num_items = len(prt_items)
            if num_items <= self.study_items_max:
                return "\n".join([", ".join(str(e) for e in sublist) for sublist in prt_multiline])
            else:
                short_list = prt_items[:self.study_items_max]
                short_mult = [short_list[i:i+N] for i in range(0, len(short_list), N)]
                short_str = "\n".join([", ".join(str(e) for e in sublist) for sublist in short_mult])
                return "".join(["{N} genes; ".format(N=num_items), short_str, "..."])

    def __get_genestr(self, itemid):
        """Given a geneid, return the string geneid or a gene symbol."""
        if self.id2symbol is not None:
            symbol = self.id2symbol.get(itemid, None)
            if symbol is not None:
                return symbol
        if isinstance(itemid, int):
            return str(itemid)
        return itemid

# Copyright (C) 2016, DV Klopfenstein, H Tang, All rights reserved.
