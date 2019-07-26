"""GO Grouper Plotting objects."""

import sys
import os
import collections as cx
from goatools.gosubdag.plot.gosubdag_plot import GoSubDagPlot
from goatools.gosubdag.gosubdag import GoSubDag
from goatools.gosubdag.go_tasks import get_go2parents_go2obj
from goatools.grouper.colors import GrouperColors

__copyright__ = "Copyright (C) 2016-2018, DV Klopfenstein, H Tang, All rights reserved."
__author__ = "DV Klopfenstein"


# pylint: disable=too-many-instance-attributes
class PltGroupedGosArgs(object):
    """Arguments for plot_groups_all."""

    keys_plt = set(['go2color', 'go2name', 'go2bordercolor', 'max_gos', 'go2txt'])
    keys_loc = set(['max_gos', 'upper_trigger', 'max_upper', 'fout_dir', 'hdrgo_dflt_color'])
    hdrgo_dflt_color = '#ffedcc' # very light pale orange

    def __init__(self, grprobj, **kws):
        # Keyword args for GO group plots
        self.kws = kws.copy()
        _objcolors = GrouperColors(grprobj)
        self.go2color = kws.get('go2color', _objcolors.get_go2color_users())
        self.go2name = kws.get('go2name', None)
        self.go2bordercolor = self._init_go2bordercolor(_objcolors, **kws)
        # Location
        self.max_gos = kws.get('max_gos', 100)
        self.upper_trigger = kws.get('upper_trigger', 50)
        self.max_upper = kws.get('max_upper', 50)
        self.fout_dir = kws.get('fout_dir', None)
        # png
        self.plt_ext = kws.get('plt_ext', "png")

    def get_go2color_inst(self, hdrgo):
        """Get a copy of go2color with GO group header colored."""
        go2color = self.go2color.copy()
        go2color[hdrgo] = self.hdrgo_dflt_color
        return go2color

    def get_kws_plt(self):
        """Get keyword args for GoSubDagPlot from self unless they are None."""
        kws_plt = {}
        for key_plt in self.keys_plt:
            key_val = getattr(self, key_plt, None)
            if key_val is not None:
                kws_plt[key_plt] = key_val
            elif key_plt in self.kws:
                kws_plt[key_plt] = self.kws[key_plt]
        return kws_plt

    @staticmethod
    def _init_go2bordercolor(objcolors, **kws):
        """Initialize go2bordercolor with default to make hdrgos bright blue."""
        go2bordercolor_ret = objcolors.get_bordercolor()
        if 'go2bordercolor' not in kws:
            return go2bordercolor_ret
        go2bordercolor_usr = kws['go2bordercolor']
        goids = set(go2bordercolor_ret).intersection(go2bordercolor_usr)
        for goid in goids:
            go2bordercolor_usr[goid] = go2bordercolor_ret[goid]
        return go2bordercolor_usr


# --------------------------------------------------------------------------------------
class PltGroupedGos(object):
    """Plots grouped user GO IDs."""

    ntpltgo = cx.namedtuple("NtPltGo", "hdrgo gosubdag tot_usrgos parentcnt desc")

    def __init__(self, grprplt):
        self.grprplt = grprplt  # GrouperPlot
        self.usrgos = grprplt.grprobj.usrgos  # set(['GO:NNNNNNN', ...
        self.gosubdag = grprplt.grprobj.gosubdag  # GoSubDag
        self.grprobj = grprplt.grprobj  # Grouper

    def plot_sections(self, fout_dir=".", **kws_pltargs):
        """Plot GO DAGs for all sections (not Misc.)."""
        hdrgos = self.grprobj.hdrobj.get_section_hdrgos()
        pltargs = PltGroupedGosArgs(self.grprobj, fout_dir=fout_dir, **kws_pltargs)
        return self._plot_groups_hdrgos(hdrgos, pltargs)

    def plot_groups_all(self, fout_dir=".", **kws_pltargs): # GrouperUserGos
        """Plot GO DAGs for all groups of user GOs."""
        hdrgos = self.grprobj.get_hdrgos()
        pltargs = PltGroupedGosArgs(self.grprobj, fout_dir=fout_dir, **kws_pltargs)
        return self._plot_groups_hdrgos(hdrgos, pltargs)

    def plot_groups_unplaced(self, fout_dir=".", **kws_pltargs):
        """Plot GO DAGs for groups of user GOs which are not in a section."""
        hdrgos = self.grprobj.get_hdrgos_unplaced()
        pltargs = PltGroupedGosArgs(self.grprobj, fout_dir=fout_dir, **kws_pltargs)
        return self._plot_groups_hdrgos(hdrgos, pltargs)

    def get_pltdotstrs(self, **kws): # GrouperUserGos
        """Plot GO DAGs for all groups of user GOs."""
        hdrgos = self.grprobj.get_hdrgos()
        return self._get_pltdotstrs(hdrgos, **kws)

    def _get_plt_data(self, hdrgos_usr):
        """Given User GO IDs, return their GO headers and other GO info."""
        hdrgo2usrgos = self.grprobj.get_hdrgo2usrgos(hdrgos_usr)
        usrgos_actual = set([u for us in hdrgo2usrgos.values() for u in us])
        go2obj = self.gosubdag.get_go2obj(usrgos_actual.union(hdrgo2usrgos.keys()))
        return hdrgo2usrgos, go2obj

    def _get_pltdotstrs(self, hdrgos_usr, **kws):
        """Plot GO DAGs for each group found under a specfied header GO."""
        import datetime
        import timeit
        dotstrs_all = []
        tic = timeit.default_timer()
        # Loop through GO groups. Each group of GOs is formed under a single "header GO"
        hdrgo2usrgos, go2obj = self._get_plt_data(hdrgos_usr)
        # get dot strings with _get_dotstrs_curs
        for hdrgo, usrgos in hdrgo2usrgos.items():
            dotstrs_cur = self._get_dotgraphs(
                hdrgo, usrgos,
                pltargs=PltGroupedGosArgs(self.grprobj, **kws),
                # TBD: get_go2parents_go2obj: Add user-specified relationships
                go2parentids=get_go2parents_go2obj(go2obj))
            dotstrs_all.extend(dotstrs_cur)
        sys.stdout.write("\nElapsed HMS: {HMS} to write ".format(
            HMS=str(datetime.timedelta(seconds=(timeit.default_timer()-tic)))))
        sys.stdout.write("{P:5,} GO DAG plots for {H:>5,} GO grouping headers\n".format(
            H=len(hdrgo2usrgos), P=len(dotstrs_all)))
        return sorted(set(dotstrs_all))

    def _plot_groups_hdrgos(self, hdrgos_usr, pltargs):
        """Plot GO DAGs for each group found under a specfied header GO."""
        import datetime
        import timeit
        pngs = []
        tic = timeit.default_timer()
        # Loop through GO groups. Each group of GOs is formed under a single "header GO"
        hdrgo2usrgos, go2obj = self._get_plt_data(hdrgos_usr)
        for hdrgo, usrgos in hdrgo2usrgos.items():
            pngs.extend(self._plot_go_group(
                hdrgo, usrgos,
                pltargs=pltargs,
                # TBD: get_go2parents_go2obj: Add user-specified relationships
                go2parentids=get_go2parents_go2obj(go2obj)))
        sys.stdout.write("\nElapsed HMS: {HMS} to write ".format(
            HMS=str(datetime.timedelta(seconds=(timeit.default_timer()-tic)))))
        sys.stdout.write("{P:5,} GO DAG plots for {H:>5,} GO grouping headers\n".format(
            H=len(hdrgo2usrgos), P=len(pngs)))
        return sorted(set(pngs))

    def _plot_go_group(self, hdrgo, usrgos, pltargs, go2parentids):
        """Plot an exploratory GO DAG for a single Group of user GOs."""
        gosubdagplotnts = self._get_gosubdagplotnts(hdrgo, usrgos, pltargs, go2parentids)
        # Create pngs and return png names
        pngs = [obj.wrplt(pltargs.fout_dir, pltargs.plt_ext) for obj in gosubdagplotnts]
        return pngs

    def _get_dotgraphs(self, hdrgo, usrgos, pltargs, go2parentids):
        """Get a GO DAG in a dot-language string for a single Group of user GOs."""
        gosubdagplotnts = self._get_gosubdagplotnts(hdrgo, usrgos, pltargs, go2parentids)
        # Create DAG graphs as dot language strings. Loop through GoSubDagPlotNt list
        dotstrs = [obj.get_dotstr() for obj in gosubdagplotnts]
        return dotstrs

    def _get_gosubdagplotnts(self, hdrgo, usrgos, pltargs, go2parentids):
        """Get list of GoSubDagPlotNt for plotting an exploratory GODAG for 1 Group of user GOs."""
        dotgraphs = []
        go2color = pltargs.get_go2color_inst(hdrgo)
        # namedtuple fields: hdrgo gosubdag tot_usrgos parentcnt desc
        ntpltgo0 = self._get_pltdag_ancestors(hdrgo, usrgos, desc="")
        ntpltgo1 = self._get_pltdag_path_hdr(hdrgo, usrgos, desc="pruned")
        num_go0 = len(ntpltgo0.gosubdag.go2obj)
        num_go1 = len(ntpltgo1.gosubdag.go2obj)
        title = "{GO} {NAME}; {N} GO sources".format(
            GO=hdrgo,
            NAME=self.gosubdag.go2obj[hdrgo].name,
            N=len(ntpltgo0.gosubdag.go_sources))
        # print("PltGroupedGos::_get_gosubdagplotnts TITLE", title)
        if num_go0 < pltargs.max_gos:
            # print("PltGroupedGos::_get_gosubdagplotnts ntpltgo0 ALWAYS IF NOT TOO BIG")
            dotgraphs.append(self._get_gosubdagplotnt(ntpltgo0, title, go2color, pltargs))
        # PLOT A: Plot the entire GO ID group under GO header, hdrgo, if not too big
        if num_go0 < pltargs.max_gos and \
            ntpltgo0.tot_usrgos == ntpltgo1.tot_usrgos:
            # print("PltGroupedGos::_get_gosubdagplotnts ntpltgo0 AAAAAAAAAAAAAAAAA")
            dotgraphs.append(self._get_gosubdagplotnt(ntpltgo0, title, go2color, pltargs))
        # PLOT B: Plot only the GO ID group passing thru the GO header, hdrgo, if not too big
        elif num_go1 < pltargs.max_gos and ntpltgo0.tot_usrgos != ntpltgo1.tot_usrgos:
            # print("PltGroupedGos::_get_gosubdagplotnts ntpltgo1(pruned) BBBBBBBBBBBBBBBBB")
            dotgraphs.append(self._get_gosubdagplotnt(ntpltgo1, title, go2color, pltargs))
        # PLOT C: If DAG is large, just print upper portion
        # PLOT D: If DAG is large, just print upper portion passing through the GO header
        elif num_go1 >= pltargs.upper_trigger:
            # print("PltGroupedGos::_get_gosubdagplotnts upper_pruned CCCCCCCCCCCCCCCCC")
            gos_upper = self._get_gos_upper(ntpltgo1, pltargs.max_upper, go2parentids)
            #ntpltgo2 = self._get_pltdag_ancestors(hdrgo, gos_upper, "{BASE}_upper.png")
            ntpltgo3 = self._get_pltdag_path_hdr(hdrgo, gos_upper, "upper_pruned")
            # Middle GO terms chosen to be start points will be green unless reset back
            for goid in gos_upper:
                if goid not in go2color:
                    go2color[goid] = 'white'
            dotgraphs.append(self._get_gosubdagplotnt(ntpltgo3, title, go2color, pltargs))
        else:
            # print("PltGroupedGos::_get_gosubdagplotnts EEEEEEEEEEEEEEEEE")
            self._no_ntplt(ntpltgo0)
        return dotgraphs

    def _get_gos_upper(self, ntpltgo1, max_upper, go2parentids):
        """Plot a GO DAG for the upper portion of a single Group of user GOs."""
        # Get GO IDs which are in the hdrgo path
        goids_possible = ntpltgo1.gosubdag.go2obj.keys()
        # Get upper GO IDs which have the most descendants
        return self._get_gosrcs_upper(goids_possible, max_upper, go2parentids)

    def _get_gosrcs_upper(self, goids, max_upper, go2parentids):
        """Get GO IDs for the upper portion of the GO DAG."""
        gosrcs_upper = set()
        get_nt = self.gosubdag.go2nt.get
        go2nt = {g:get_nt(g) for g in goids}
        # Sort by descending order of descendant counts to find potential new hdrgos
        go_nt = sorted(go2nt.items(), key=lambda t: -1*t[1].dcnt)
        goids_upper = set()
        for goid, _ in go_nt: # Loop through GO ID, GO nt
            goids_upper.add(goid)
            if goid in go2parentids:
                goids_upper |= go2parentids[goid]
            #print "{} {:3} {}".format(goid, len(goids_upper), gont.GO_name)
            if len(goids_upper) < max_upper:
                gosrcs_upper.add(goid)
            else:
                break
        return gosrcs_upper

    def _get_gosubdagplotnt(self, ntplt, title, go2color, pltargs):
        """Return GoSubDagPlotNt, which contains both a GoSubDagPlot object and ntobj."""
        kws_plt = pltargs.get_kws_plt()
        kws_plt['id'] = '"{ID}"'.format(ID=ntplt.hdrgo)
        kws_plt['title'] = "{TITLE} of {M} user GOs".format(TITLE=title, M=ntplt.tot_usrgos)
        kws_plt['go2color'] = go2color
        kws_plt['go2bordercolor'] = pltargs.go2bordercolor
        if ntplt.parentcnt:
            kws_plt["parentcnt"] = True
        gosubdagplot = GoSubDagPlot(ntplt.gosubdag, **kws_plt)
        return GoSubDagPlotNt(self.grprobj, gosubdagplot, ntplt)

    def _no_ntplt(self, ntplt):
        """Print a message about the GO DAG Plot we are NOT plotting."""
        sys.stdout.write("  {GO_USR:>6,} usr {GO_ALL:>6,} GOs  DID NOT WRITE: {B} {D}\n".format(
            B=self.grprobj.get_fout_base(ntplt.hdrgo),
            D=ntplt.desc,
            GO_USR=len(ntplt.gosubdag.go_sources),
            GO_ALL=len(ntplt.gosubdag.go2obj)))

    def _get_pltdag_ancestors(self, hdrgo, usrgos, desc=""):
        """Get GoSubDag containing hdrgo and all usrgos and their ancestors."""
        go_srcs = usrgos.union([hdrgo])
        gosubdag = GoSubDag(go_srcs,
                            self.gosubdag.get_go2obj(go_srcs),
                            relationships=self.gosubdag.relationships,
                            rcntobj=self.gosubdag.rcntobj,
                            go2nt=self.gosubdag.go2nt)
        tot_usrgos = len(set(gosubdag.go2obj.keys()).intersection(self.usrgos))
        return self.ntpltgo(
            hdrgo=hdrgo,
            gosubdag=gosubdag,
            tot_usrgos=tot_usrgos,
            parentcnt=False,
            desc=desc)

    def _get_pltdag_path_hdr(self, hdrgo, usrgos, desc="pruned"):
        """Get GoSubDag with paths from usrgos through hdrgo."""
        go_sources = usrgos.union([hdrgo])
        gosubdag = GoSubDag(go_sources,
                            self.gosubdag.get_go2obj(go_sources),
                            relationships=self.gosubdag.relationships,
                            rcntobj=self.gosubdag.rcntobj,
                            go2nt=self.gosubdag.go2nt,
                            dst_srcs_list=[(hdrgo, usrgos), (None, set([hdrgo]))])
        tot_usrgos = len(set(gosubdag.go2obj.keys()).intersection(self.usrgos))
        return self.ntpltgo(
            hdrgo=hdrgo,
            gosubdag=gosubdag,
            tot_usrgos=tot_usrgos,
            parentcnt=True,
            desc=desc)

# --------------------------------------------------------------------------------------
class GoSubDagPlotNt(object):
    """Contains Grouper object, GoSubDagPlot object, and namedtuple."""

    def __init__(self, grprobj, gosubdagplot, ntplt):
        self.grprobj = grprobj
        self.gosubdagplot = gosubdagplot # GoSubDagPlot
        # namedtuple fields: hdrgo gosubdag tot_usrgos parentcnt desc
        self.ntplt = ntplt

    def wrplt(self, fout_dir, plt_ext="png"):
        """Write png containing plot of GoSubDag."""
        # Ex basename
        basename = self.grprobj.get_fout_base(self.ntplt.hdrgo)
        plt_pat = self.get_pltpat(plt_ext)
        fout_basename = plt_pat.format(BASE=basename)
        fout_plt = os.path.join(fout_dir, fout_basename)
        self.gosubdagplot.plt_dag(fout_plt) # Create Plot
        return fout_plt

    def get_dotstr(self):
        """Return a string containing DAG graph in Grpahviz's dot language."""
        dotobj = self.gosubdagplot.get_pydot_graph() # pydot.Dot
        dotstr = dotobj.create_dot()
        return dotstr

    def get_pltpat(self, plt_ext="svg"):
        """Return png pattern: {BASE}.png {BASE}_pruned.png {BASE}_upper_pruned.png"""
        if self.ntplt.desc == "":
            return ".".join(["{BASE}", plt_ext])
        return "".join(["{BASE}_", self.ntplt.desc, ".", plt_ext])


# Copyright (C) 2016-2018, DV Klopfenstein, H Tang, All rights reserved.
