"""Given user GO ids and parent terms, group user GO ids under one parent term.

   Given a group of GO ids with one or more higher-level grouping terms, group
   each user GO id under the most descriptive parent GO term.

   Each GO id may have more than one parent.  One of the parent(s) is chosen
   to best represent the user GO id's function. The choice of parent is made by
   regarding how close the parent GO id is to the bottom of its hierarchy.

   The estimation of how close a GO term is to "the bottom" of its GO hierarchy
   is estimated using the number of total Go term descendent counts below
   that term.
"""

from goatools.gosubdag.plot.gosubdag_plot import GoSubDagPlot
from goatools.gosubdag.gosubdag import GoSubDag
from goatools.grouper.plotobj import PltGroupedGos, PltGroupedGosArgs
from goatools.grouper.colors import GrouperColors
from goatools.grouper.grprobj import Grouper

__copyright__ = "Copyright (C) 2016-2018, DV Klopfenstein, H Tang, All rights reserved."
__author__ = "DV Klopfenstein"


class GrouperPlot(object):
    """Groups the user GO ids under other GO IDs acting as headers for the GO groups."""

    def __init__(self, grprobj):
        self.grprobj = grprobj

    def plot_sections(self, fout_dir=".", **kws_usr):
        """Plot groups of GOs which have been placed in sections."""
        kws_plt, _ = self._get_kws_plt(None, **kws_usr)
        PltGroupedGos(self).plot_sections(fout_dir, **kws_plt)

    def plot_groups_all(self, fout_dir=".", **kws_pltargs): # Grouper
        """Plot each GO header group in Grouper."""
        # kws: go2color max_gos upper_trigger max_upper
        return PltGroupedGos(self).plot_groups_all(fout_dir, **kws_pltargs)

    def get_pltdotstrs(self, **kws_usr):
        """Plot each GO header group in Grouper."""
        # kws: go2color max_gos upper_trigger max_upper
        return PltGroupedGos(self).get_pltdotstrs(**kws_usr)

    def get_pltdotstr(self, **kws_usr):
        """Plot one GO header group in Grouper."""
        dotstrs = self.get_pltdotstrs(**kws_usr)
        assert len(dotstrs) == 1
        return dotstrs[0]

    def plot_groups_unplaced(self, fout_dir=".", **kws_usr):
        """Plot each GO group."""
        # kws: go2color max_gos upper_trigger max_upper
        plotobj = PltGroupedGos(self)
        return plotobj.plot_groups_unplaced(fout_dir, **kws_usr)

    def get_gosubdagplot(self, goids=None, **kws_usr):
        """Plot GO IDs."""
        if goids is None:
            goids = self.grprobj.usrgos
        kws_plt, kws_dag = self._get_kws_plt(goids, **kws_usr)
        gosubdag = GoSubDag(
            goids,
            self.grprobj.gosubdag.get_go2obj(goids),
            self.grprobj.gosubdag.relationships,
            rcntobj=self.grprobj.gosubdag.rcntobj,
            go2nt=self.grprobj.gosubdag.go2nt,
            **kws_dag)
        return GoSubDagPlot(gosubdag, **kws_plt)

    def plot_gos(self, fout_img, goids=None, **kws_usr):
        """Plot GO IDs."""
        gosubdagplot = self.get_gosubdagplot(goids, **kws_usr)  # GoSubDagPlot
        gosubdagplot.plt_dag(fout_img)

    def _get_kws_plt(self, usrgos, **kws_usr):
        """Add go2color and go2bordercolor relevant to this grouping into plot."""
        kws_plt = kws_usr.copy()
        kws_dag = {}
        hdrgo = kws_plt.get('hdrgo', None)
        objcolor = GrouperColors(self.grprobj)
        # GO term colors
        if 'go2color' not in kws_usr:
            kws_plt['go2color'] = objcolor.get_go2color_users()
        elif hdrgo is not None:
            go2color = kws_plt.get('go2color').copy()
            go2color[hdrgo] = PltGroupedGosArgs.hdrgo_dflt_color
            kws_plt['go2color'] = go2color
        # GO term border colors
        if 'go2bordercolor' not in kws_usr:
            kws_plt['go2bordercolor'] = objcolor.get_bordercolor()
        prune = kws_usr.get('prune', None)
        if prune is True and hdrgo is not None:
            kws_dag['dst_srcs_list'] = [(hdrgo, usrgos), (None, set([hdrgo]))]
            kws_plt['parentcnt'] = True
        elif prune:
            kws_dag['dst_srcs_list'] = prune
            kws_plt['parentcnt'] = True
        # Group text
        kws_plt['go2txt'] = self.get_go2txt(self.grprobj,
                                            kws_plt.get('go2color'), kws_plt.get('go2bordercolor'))
        return kws_plt, kws_dag

    @staticmethod
    def get_go2txt(grprobj_cur, grp_go2color, grp_go2bordercolor):
        """Adds section text in all GO terms if not Misc. Adds Misc in terms of interest."""
        goids_main = set(o.id for o in grprobj_cur.gosubdag.go2obj.values())
        hdrobj = grprobj_cur.hdrobj
        grprobj_all = Grouper("all",
                              grprobj_cur.usrgos.union(goids_main), hdrobj, grprobj_cur.gosubdag)
        # Adds section text to all GO terms in plot (misses middle GO terms)
        _secdflt = hdrobj.secdflt
        _hilight = set(grp_go2color.keys()).union(grp_go2bordercolor)
        ret_go2txt = {}
        # Keep sections text only if GO header, GO user, or not Misc.
        if hdrobj.sections:
            for goid, txt in grprobj_all.get_go2sectiontxt().items():
                if txt == 'broad':
                    continue
                if txt != _secdflt or goid in _hilight:
                    ret_go2txt[goid] = txt
        return ret_go2txt

    def plot_grouped_gos(self, fout_img=None, exclude_hdrs=None, **kws_usr):
        """One Plot containing all user GOs (yellow or green) and header GO IDs(green or purple)."""
        # kws_plt -> go2color go2bordercolor
        kws_plt, kws_dag = self._get_kws_plt(self.grprobj.usrgos, **kws_usr)
        pltgosusr = self.grprobj.usrgos
        if exclude_hdrs is not None:
            pltgosusr = pltgosusr.difference(self.grprobj.get_usrgos_g_hdrgos(exclude_hdrs))
        if fout_img is None:
            fout_img = "{GRP_NAME}.png".format(GRP_NAME=self.grprobj.grpname)
        # Split one plot into potentially three (BP, MF, CC) if png filename contains '{NS}'
        if '{NS}' in fout_img:
            go2nt = self.grprobj.gosubdag.get_go2nt(pltgosusr)
            for namespace in ['BP', 'MF', 'CC']:
                pltgos_ns = [go for go in pltgosusr if go2nt[go].NS == namespace]
                if pltgos_ns:
                    png = fout_img.format(NS=namespace)
                    self._plot_grouped_gos(png, pltgos_ns, kws_plt, kws_dag)
        # Plot all user GO IDs into a single plot, regardless of their namespace
        else:
            self._plot_grouped_gos(fout_img, pltgosusr, kws_plt, kws_dag)

    def _plot_grouped_gos(self, fout_img, pltgosusr, kws_plt, kws_dag):
        gosubdag_plt = GoSubDag(
            pltgosusr,
            self.grprobj.gosubdag.get_go2obj(pltgosusr),
            self.grprobj.gosubdag.relationships,
            rcntobj=self.grprobj.gosubdag.rcntobj,
            go2nt=self.grprobj.gosubdag.go2nt,
            **kws_dag)
        godagplot = GoSubDagPlot(gosubdag_plt, **kws_plt)
        godagplot.plt_dag(fout_img)

# Copyright (C) 2016-2018, DV Klopfenstein, H Tang, All rights reserved.
