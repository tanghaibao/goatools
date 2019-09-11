"""Plot a GoSubDag.

   GO Terms in a plot contain text. The first two lines appear by default:

        GO:0015618 L07 D13 p2 c2 d3
        potassium-transporting ATPase activity

   GO:0015618 => GO ID

   First line:

   LNN        => The "level" of the go term.
                 The length of the shortest path(s) from the top to the current term.

   DNN        => The "depth" of the go term.
                 The length of the longest path(s) from the top to the current term.

   pN         => Optional (parentcnt arg): number of immediate parent terms if the number
                 of parents plotted is different than the number of parents in the obo.
                 The default is to plot all higher-level parent terms from the
                 source GO IDs to best show the hierarchy.
                 But with some plots, plotting all parent GO terms can result in
                 a GO DAG which is too large to be clearly readable. The user can then
                 force some parent GO terms to not be plotting using the init_dag=path keyword.
                 If some parent terms are not plotted, using the "parentcnt" option
                 can help the user know where the plot was cut for readability.

   cN         => Optional (childcnt arg): number of immediate child terms.
                 Child hierarchy is not traversed. The default is to not plot all lower-level
                 child terms to prevent the GO plots from being massive and unreadable.
                 So knowing the total number of immediate child terms present (with
                 most not on the plot) can give the user a better sense of the qualities
                 of their plot.

   dN         => Optional (rcntobj arg): total number of all levels of child terms.
                 Child hierarchy is traversed to the bottom or the leaf-level of the graph.
                 Knowing the total number of all descendants succinctly gives the user
                 a sense of how close a GO term is to the bottom of the graph.
                 "Descendants Count" is used a proxy for understanding if the
                 GO term is a "broad" or "specific". If the GO term broadly
                 describes a biological process, it most always has hundreds or thousands
                 of total child terms. If the GO term specifically describes
                 a biological process, it often has tens or less of total child terms.

"""

from __future__ import print_function

__copyright__ = "Copyright (C) 2016-2020, DV Klopfenstein, H Tang, All rights reserved."
__author__ = "DV Klopfenstein"

import sys
import os
import pydot
from goatools.gosubdag.go_edges import get_edgesobj
from goatools.gosubdag.plot.go_node import GoNodeOpts
from goatools.gosubdag.plot.go_node import GoNode
from goatools.gosubdag.plot.go_edge import GoEdgeOpts
from goatools.gosubdag.plot.go_edge import GoEdge
from goatools.gosubdag.plot.go2color import Go2Color
from goatools.gosubdag.plot.goea_results import GoeaResults
from goatools.gosubdag.utils import get_kwargs


class GoSubDagPlot:
    """Plot a graph contained in an object of type GoSubDag ."""

    # http://www.graphviz.org/doc/info/colors.html
    # pylint: disable=bad-whitespace
    rel2edgekws = {
        'is_a':                 {'color':'black',       'style':'solid'},
        'part_of':              {'color':'#ff5b00',     'style':'dashed'},  # xkcd bright orange
        'regulates':            {'color':'purple3',     'style':'dashed'},
        'positively_regulates': {'color':'magenta',     'style':'dashed'},
        'negatively_regulates': {'color':'#41fdfe',     'style':'dashed'},  # xkcd bright cyan
        'occurs_in':            {'color':'aquamarine4', 'style':'dashed'},
        'capable_of':           {'color':'dodgerblue',  'style':'dashed'},
        'capable_of_part_of':   {'color':'darkorange',  'style':'dashed'},
    }

    # https://graphviz.org/doc/info/attrs.html
    exp_keys = {
        'dag': {'title', 'id', 'dpi'},  # pydot.Dot kwargs
        # goobj2fncname parentcnt shorten mark_alt_id childcnt prt_pcnt ...
        'node_go': GoNodeOpts.exp_keys.union(GoNodeOpts.exp_elems),
        'edge_go': {'edge2txt'},
        # id2symbol study_items items_p_line pval_name
        'goea': GoeaResults.kws_set,
    }

    dflts = {'dpi':150}

    def __init__(self, gosubdag, **kwu):
        assert gosubdag, "**FATAL: MISSING SUBSET GODag"
        self.gosubdag = gosubdag
        # kwu: id, title, dpi, go2txt
        self.kws = self._init_kws(**kwu)
        # kwu: log parentcnt
        self.edgesobj = get_edgesobj(gosubdag, **kwu)
        # pylint: disable=line-too-long
        # kwu: go2color go2bordercolor dflt_bordercolor
        _node_opt = kwu['GoNodeOpts'] if 'GoNodeOpts' in kwu else self._init_gonodeopts(**kwu)
        _edge_opt = kwu['GoEdgeOpts'] if 'GoEdgeOpts' in kwu else self._init_goedgeopts()
        _objcolor = kwu['Go2Color'] if 'Go2Color' in kwu else self._init_objcolor(_node_opt, **kwu)
        self.pydotnodego = GoNode(gosubdag, _objcolor, _node_opt)
        self.pydotedge = GoEdge(gosubdag, _edge_opt)
        self.log = kwu.get('log', sys.stdout)
        #    KWS=kws.keys(), V=kws['parentcnt'] if 'parentcnt' in kws else None)

    def _init_objcolor(self, node_opts, **kwu):
        """Return user-created Go2Color object or create one."""
        objgoea = node_opts.kws['dict'].get('objgoea', None)
        # kwu: go2color go2bordercolor dflt_bordercolor key2col
        return Go2Color(self.gosubdag, objgoea, **kwu)

    def _init_gonodeopts(self, **kws_usr):
        """Initialize a GO Node plot options object, GoNodeOpts."""
        options = GoNodeOpts(self.gosubdag, **self.kws['node_go'])
        # Add parent edge count if either is in kws: parentcnt, prt_pcnt
        if not options.kws['set'].isdisjoint(['parentcnt', 'prt_pcnt']):
            options.kws['dict']['c2ps'] = self.edgesobj.get_c2ps()
        # GoeaResults(kws['goea_results'], **self.kws['goea']) if 'goea_results' in kws else None
        if 'goea_results' in kws_usr:
            objgoea = GoeaResults(kws_usr['goea_results'], **self.kws['goea'])
            options.kws['dict']['objgoea'] = objgoea
        return options

    def _init_goedgeopts(self):
        """Initialize a GO Edge plot options object, GoEdgeOpts."""
        options = GoEdgeOpts(self.gosubdag, **self.kws['edge_go'])
        return options

    def prt_goids(self, prt):
        """Print all GO IDs in the plot, plus their color."""
        fmt = self.gosubdag.prt_attr['fmta']
        nts = sorted(self.gosubdag.go2nt.values(), key=lambda nt: [nt.NS, nt.depth, nt.alt])
        _get_color = self.pydotnodego.go2color.get
        for ntgo in nts:
            gostr = fmt.format(**ntgo._asdict())
            col = _get_color(ntgo.GO, "")
            prt.write("{COLOR:7} {GO}\n".format(COLOR=col, GO=gostr))

    def _init_kws(self, **kws_usr):
        """Return a dict containing user-specified plotting options."""
        kws_self = {}
        user_keys = set(kws_usr)
        for objname, expset in self.exp_keys.items():
            usrkeys_curr = user_keys.intersection(expset)
            kws_self[objname] = get_kwargs(kws_usr, usrkeys_curr, usrkeys_curr)
        if 'title' in kws_usr:
            kws_self['dag']['label'] = kws_usr['title']
            kws_self['dag']['labelloc'] = 't'
        dpi = str(kws_self['dag'].get('dpi', self.dflts['dpi']))
        kws_self['dag']['dpi'] = dpi
        return kws_self

    def plt_dag(self, fout_img, engine="pydot"):
        """Plot using pydot, graphviz, or GML."""
        if engine == "pydot":
            self._plt_pydot(fout_img)
        else:
            raise RuntimeError("ENGINE NOT IMPLEMENTED({E})".format(E=engine))

    # ----------------------------------------------------------------------------------
    # pydot
    def _plt_pydot(self, fout_img):
        """Plot using the pydot graphics engine."""
        dag = self.get_pydot_graph()
        self.wr_pydot_dag(fout_img, dag)

    def wr_pydot_dag(self, fout_img, dag):
        """Plot using the pydot graphics engine."""
        img_fmt = os.path.splitext(fout_img)[1][1:]
        dag.write(fout_img, format=img_fmt)
        self.log.write("  {GO_USR:>3} usr {GO_ALL:>3} GOs  WROTE: {F}\n".format(
            F=fout_img,
            GO_USR=len(self.gosubdag.go_sources),
            GO_ALL=len(dag.obj_dict['nodes'])))

    def get_pydot_graph(self):
        """Given a DAG, return a pydot digraph object."""
        rel = "is_a"
        # Initialize empty dag
        dag = pydot.Dot(graph_type='digraph', **self.kws['dag'])
        # Initialize nodes
        go2node = self._get_go2pydotnode()
        # Add nodes to graph
        for node in go2node.values():
            dag.add_node(node)
        # Add edges to graph
        rel2edgekws = self.rel2edgekws
        self.edgesobj.chk_edges()
        edgekws = rel2edgekws.get(rel)
        _add_edges = self.pydotedge.add_edges
        _add_edges(self.edgesobj.edges, go2node, dag, **edgekws)
        for reltype, edges_list in self.edgesobj.edges_rel.items():
            edgekws = rel2edgekws.get(reltype)
            _add_edges(edges_list, go2node, dag, **edgekws)
        return dag

    @staticmethod
    def _prt_edge(dag_edge, attr):
        """Print edge attribute"""
        # sequence parent_graph points attributes type parent_edge_list
        print("Edge {ATTR}: {VAL}".format(ATTR=attr, VAL=dag_edge.obj_dict[attr]))

    def _get_go2pydotnode(self):
        """Create pydot Nodes."""
        go2node = {}
        go2obj = self.gosubdag.go2obj
        get_node = self.pydotnodego.get_node
        for goid in self.get_goids_plt():
            goobj = go2obj[goid]
            node = get_node(goid, goobj)
            go2node[goid] = node
            if goid != goobj.id: # goid is an alias for goobj
                go2node[goobj.id] = node
        return go2node

    def get_goids_plt(self):
        """Get GO IDs to be plotted, given a GoSubDag."""
        return self.edgesobj.get_all_edge_nodes()


# Copyright (C) 2016-2020, DV Klopfenstein, H Tang, All rights reserved.
