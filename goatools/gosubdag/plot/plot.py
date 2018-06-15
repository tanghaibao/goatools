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

__copyright__ = "Copyright (C) 2016-2017, DV Klopfenstein, H Tang, All rights reserved."
__author__ = "DV Klopfenstein"

import collections as cx
from goatools.gosubdag.gosubdag import GoSubDag
from goatools.gosubdag.plot.gosubdag_plot import GoSubDagPlot


def plt_goids(gosubdag, fout_img, goids, **kws_plt):
    """Plot GO IDs in a DAG (Directed Acyclic Graph)."""
    gosubdag_plt = GoSubDag(goids, gosubdag.go2obj, rcntobj=gosubdag.rcntobj, **kws_plt)
    godagplot = GoSubDagPlot(gosubdag_plt, **kws_plt)
    godagplot.plt_dag(fout_img)
    return godagplot

def plot_gos(fout_img, goids, go2obj, **kws):
    """Given GO ids and the obo_dag, create a plot of paths from GO ids."""
    gosubdag = GoSubDag(goids, go2obj, rcntobj=True)
    godagplot = GoSubDagPlot(gosubdag, **kws)
    godagplot.plt_dag(fout_img)

def plot_results(fout_img, goea_results, **kws):
    """Given a list of GOEA results, plot result GOs up to top."""
    if "{NS}" not in fout_img:
        plt_goea_results(fout_img, goea_results, **kws)
    else:
        # Plot separately by NS: BP, MF, CC
        ns2goea_results = cx.defaultdict(list)
        for rec in goea_results:
            ns2goea_results[rec.NS].append(rec)
        for ns_name, ns_res in ns2goea_results.items():
            fout = fout_img.format(NS=ns_name)
            plt_goea_results(fout, ns_res, **kws)

def plt_goea_results(fout_img, goea_results, **kws):
    """Plot a single page."""
    go_sources = [rec.GO for rec in goea_results]
    go2obj = {rec.GO:rec.goterm for rec in goea_results}
    gosubdag = GoSubDag(go_sources, go2obj, rcntobj=True)
    godagplot = GoSubDagPlot(gosubdag, goea_results=goea_results, **kws)
    godagplot.plt_dag(fout_img)


# Copyright (C) 2016-2017, DV Klopfenstein, H Tang, All rights reserved.
