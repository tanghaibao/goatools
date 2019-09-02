"""Command-line script to create GO term diagrams

Usage:
  go_plot.py [GO ...] [options]
  go_plot.py [GO ...] [--obo=<file.obo>] [--outfile=<file.png>] [--title=<title>]
             [--go_file=<file.txt>]
             [--relationship]
             [--relationships=<part_of>]
             [--sections=<sections.txt>]
             [--gpad=<file.gpad>]
             [--gaf=<file.gaf>]
             [--gene2go=<gene2go>] [--taxid=<Taxonomy_number>]
             [--id2gos=<file.txt>]
             [--shorten]
             [--parentcnt] [--childcnt] [--mark_alt_id]
             [--go_aliases=<go_aliases.txt>]
             [--draw-children]
             [--norel]
             [--go_color_file=<file.txt>]
  go_plot.py [GO ...] [--obo=<file.obo>] [-o <file.png>] [-t <title>]
             [--shorten] [-p] [-c]
  go_plot.py [GO ...] [-o <file.png>] [--draw-children]
  go_plot.py [GO ...] [-o <file.png>] [--draw-children] [--shorten]
  go_plot.py [--obo=<file.obo>]
  go_plot.py [--obo=<file.obo>] [--outfile=<file.png>]
  go_plot.py [GO ...]
  go_plot.py [GO ...] [--outfile=<file.png>] [--title=<title>]
  go_plot.py [GO ...] [--outfile=<file.png>] [--title=<title>] [--shorten]
  go_plot.py [GO ...] [-o <file.png>] [-t <title>]
  go_plot.py [GO ...] [-o <file.png>] [-t <title>] [--parentcnt]
  go_plot.py [GO ...] [-o <file.png>] [-t <title>] [--childcnt]
  go_plot.py [GO ...] [-o <file.png>] [-t <title>] [--parentcnt] [--childcnt]
  go_plot.py [GO ...] [-o <file.png>] [-t <title>] [-p]
  go_plot.py [GO ...] [-o <file.png>] [-t <title>] [-p] [-c]

Options:
  -h --help                            show this help message and exit
  -i --go_file=<file.txt>              GO IDs in an ASCII file
  -o <file.png>, --outfile=<file.png>  Plot file name [default: go_plot.png]
  -r --relationship                    Plot all relationships
  --relationships=<part_of>            Plot user-specified relationships
  -s <sections.txt> --sections=<sections.txt>  Sections file for grouping
  -S <sections module str>             Sections file for grouping

  --gpad=<file.gpad>                   Annotations from a gpad file
  --gaf=<file.gaf>                     Annotations from a gaf file
  --id2gos=<file.txt>                  Annotations from a text file, e.g., data/association
  --gene2go=<gene2go>                  Annotations from a gene2go file downloaded from NCBI
  --taxid=<taxid_num>                  TaxID for use when reading NCBI's gene2go file

  --obo=<file.obo>                     Ontologies in obo file [default: go-basic.obo].

  -t <title>, --title=<title>          Title string to place in image
  -p --parentcnt                       Include parent count in each GO term
  -c --childcnt                        Include child count in each GO term
  --shorten                            Shorten the GO name on plots
  --mark_alt_id                        Add 'a' if GO ID is an alternate ID: GO:0007582a
  --draw-children                      Draw children. By default, they are not drawn.
  --go_aliases=<go_aliases.txt>        ASCII file containing letter alias
  --go_color_file=<file.txt>           GO color file. GO and color (eg #cafffb)

  --norel                              Don't load relationship from the GO DAG
"""

from __future__ import print_function

__copyright__ = "Copyright (C) 2016-2019, DV Klopfenstein, H Tang. All rights reserved."
__author__ = "DV Klopfenstein"


import re
import os
import sys

from goatools.obo_parser import GODag
from goatools.associations import get_tcntobj
from goatools.godag.obo_optional_attributes import OboOptionalAttrs
from goatools.godag.consts import RELATIONSHIP_SET
from goatools.godag.consts import chk_relationships
from goatools.anno.annoreader_base import AnnoReaderBase

from goatools.cli.docopt_parse import DocOptParse
from goatools.cli.gos_get import GetGOs as ClassGetGOs

from goatools.gosubdag.plot.gosubdag_plot import GoSubDagPlot
from goatools.gosubdag.plot.go2color import Go2Color
from goatools.gosubdag.gosubdag import GoSubDag
from goatools.gosubdag.go_tasks import get_go2obj_unique
from goatools.gosubdag.go_tasks import get_leaf_children
from goatools.gosubdag.rpt.wr_xlsx import read_d1_letter

from goatools.grouper.read_goids import read_sections
from goatools.grouper.grprdflts import GrouperDflts
from goatools.grouper.hdrgos import HdrgosSections
from goatools.grouper.grprobj import Grouper
from goatools.grouper.colors import GrouperColors
from goatools.grouper.grprplt import GrouperPlot


# pylint: disable=too-few-public-methods
class CliGetGOs(object):
    """Return a list of GO IDs for plotting."""

    max_gos = 200  # Maximum number of source GO IDs

    def __init__(self, go2obj):
        self.go2obj = go2obj

    def get_go_color(self, **kws):
        """Return source GO IDs and GO color, if provided."""
        ret = {'GOs':set(), 'go2color':{}}
        if 'go_color_file' in kws:
            _, go2color = ClassGetGOs.rdtxt_gos_color(kws['go_color_file'])
            self._update_ret(ret, None, go2color)
        if 'GO' in kws:
            # goids, go2color = self._goargs(ret, kws['GO'])
            goids, go2color = ClassGetGOs.get_goargs(kws['GO'], prt=sys.stdout)
            self._update_ret(ret, goids, go2color)
        if 'go_file' in kws:
            goids, go2color = ClassGetGOs.rdtxt_gos_color(kws['go_file'])
            self._update_ret(ret, goids, go2color)
        if 'draw-children' in kws:
            ret['GOs'].update(get_leaf_children(ret['GOs'], self.go2obj))
        # If there have been no GO IDs explicitly specified by the user
        if not ret['GOs']:
            # If the GO-DAG is sufficiently small, print all GO IDs
            if len(self.go2obj) < self.max_gos:
                main_gos = set(o.id for go, o in self.go2obj.items() if go != o.id)
                go_leafs = set(go for go, o in self.go2obj.items() if not o.children)
                ret['GOs'] = go_leafs.difference(main_gos)
        go2obj = {go:self.go2obj[go] for go in ret['GOs']}
        ret['GOs'] = set(get_go2obj_unique(go2obj))
        return [ret['GOs'], ret['go2color']]

    @staticmethod
    def _update_ret(ret, goids, go2color):
        """Update 'GOs' and 'go2color' in dict with goids and go2color."""
        if goids:
            ret['GOs'].update(goids)
        if go2color:
            for goid, color in go2color.items():
                ret['go2color'][goid] = color

# pylint: disable=line-too-long
class PlotCli(object):
    """Class for command-line interface for creating GO term diagrams"""

    kws_dict = set(['GO', 'outfile', 'go_file', 'sections', 'S',
                    'gpad', 'gaf', 'gene2go', 'taxid', 'id2gos',
                    'title',
                    'obo',
                    'relationships',
                    'go_color_file',
                    'go_aliases'])
    kws_set = set(['relationship',
                   'parentcnt', 'childcnt', 'mark_alt_id', 'shorten',
                   'draw-children',
                   'norel'])
    dflt_outfile = "go_plot.png"
    kws_plt = set(['parentcnt', 'childcnt', 'mark_alt_id', 'shorten'])

    def __init__(self, gosubdag=None, use_doc=True):
        _doc = __doc__ if use_doc else None
        self.objdoc = DocOptParse(_doc, self.kws_dict, self.kws_set)
        self.gosubdag = None if gosubdag is None else gosubdag

    def cli(self, kws_plt=None):
        """Command-line interface for go_draw script."""
        kws_all = self.get_docargs() if not kws_plt else kws_plt
        godag_optional_attrs = self._get_optional_attrs(kws_all)
        godag = GODag(kws_all['obo'], godag_optional_attrs)
        self.plot(godag, kws_all)

    def plot(self, godag, kws_plt):
        """Plot GO DAG subset"""
        objplt = self.get_gosubdagplot(godag, kws_plt)
        fout_img = self.get_outfile(kws_plt.get('outfile'), objplt.gosubdag.go_sources)
        objplt.prt_goids(sys.stdout)
        objplt.plt_dag(fout_img)
        #### sys.stdout.write("{N:>6} sections read\n".format(
        ####     N="NO" if sections is None else len(sections)))
        return fout_img, objplt

    def get_gosubdagplot(self, godag, kws_plt):
        """Get GoSubDagPlot"""
        # GO kws_plt: GO go_file draw-children
        goids, go2color = CliGetGOs(godag).get_go_color(**kws_plt)
        assert goids, "GO IDs NEEDED"
        #### self.gosubdag = GoSubDag(goids, godag, relationships, tcntobj=tcntobj)
        kws_dag = self._get_kwsdag(goids, godag, **kws_plt)
        relationships = self._get_relationships(kws_plt, hasattr(next(iter(godag.values())), 'relationship'))
        ## print('RRRRRRRRRRRRRRRRRRRRRRRRRR relationships', relationships)
        self.gosubdag = GoSubDag(goids, godag, relationships, **kws_dag)
        # objplt = self._plt_gogrouped(goids, go2color, **kws_plt) if 'sections' in kws_plt self._plt_gosubdag(goids, go2color, **kws_plt)
        obj = self._get_objpltg(goids, go2color, **kws_plt) if 'sections' in kws_plt else self._get_objplt(go2color, **kws_plt)
        # print('############ {N} GO IDs: relationships={Rs}'.format(N=len(obj.gosubdag.go2obj), Rs=obj.gosubdag.relationships))
        return obj

        #### if 'sections' in kws_plt:
        ####     return self._plt_gogrouped(goids, go2color, **kws_plt)
        #### else:
        ####     return self._plt_gosubdag(goids, go2color, **kws_plt)

    #### def _plt_gogrouped(self, goids, go2color_usr, **kws):
    def _get_objpltg(self, goids, go2color_usr, **kws):
        """Plot grouped GO IDs."""
        #### fout_img = self.get_outfile(kws['outfile'], goids)
        sections = read_sections(kws['sections'], exclude_ungrouped=True)
        # print ("KWWSSSSSSSS", kws)
        # kws_plt = {k:v for k, v in kws.items if k in self.kws_plt}
        grprobj_cur = self._get_grprobj(goids, sections)
        # GO: purple=hdr-only, green=hdr&usr, yellow=usr-only
        # BORDER: Black=hdr Blu=hdr&usr
        grpcolor = GrouperColors(grprobj_cur)  # get_bordercolor get_go2color_users
        grp_go2color = grpcolor.get_go2color_users()
        grp_go2bordercolor = grpcolor.get_bordercolor()
        for goid, color in go2color_usr.items():
            grp_go2color[goid] = color
        objcolor = Go2Color(self.gosubdag, objgoea=None,
                            go2color=grp_go2color, go2bordercolor=grp_go2bordercolor)
        go2txt = GrouperPlot.get_go2txt(grprobj_cur, grp_go2color, grp_go2bordercolor)
        return GoSubDagPlot(self.gosubdag, Go2Color=objcolor, go2txt=go2txt, **kws)
        #### objplt = GoSubDagPlot(self.gosubdag, Go2Color=objcolor, go2txt=go2txt, **kws)
        #### objplt.prt_goids(sys.stdout)
        #### objplt.plt_dag(fout_img)
        #### sys.stdout.write("{N:>6} sections read\n".format(
        ####     N="NO" if sections is None else len(sections)))
        #### return fout_img

    def _get_grprobj(self, goids, sections):
        """Get Grouper, given GO IDs and sections."""
        grprdflt = GrouperDflts(self.gosubdag, "goslim_generic.obo")
        hdrobj = HdrgosSections(self.gosubdag, grprdflt.hdrgos_dflt, sections)
        return Grouper("sections", goids, hdrobj, self.gosubdag)

    #### def _plt_gosubdag(self, goids, go2color, **kws):
    #### def _get_objplt(self, goids, go2color, **kws):
    def _get_objplt(self, go2color, **kws):
        """Plot GO IDs."""
        #### fout_img = self.get_outfile(kws['outfile'], goids)
        objcolor = Go2Color(self.gosubdag, objgoea=None, go2color=go2color)
        return GoSubDagPlot(self.gosubdag, Go2Color=objcolor, **kws)
        #### objplt = GoSubDagPlot(self.gosubdag, Go2Color=objcolor, **kws)
        #### objplt.prt_goids(sys.stdout)
        #### objplt.plt_dag(fout_img)
        #### return fout_img

    def _get_kwsdag(self, goids, go2obj, **kws_all):
        """Get keyword args for a GoSubDag."""
        kws_dag = {}
        # Term Counts for GO Term information score
        tcntobj = self._get_tcntobj(goids, go2obj, **kws_all)  # TermCounts or None
        if tcntobj is not None:
            kws_dag['tcntobj'] = tcntobj
        # GO letters specified by the user
        if 'go_aliases' in kws_all:
            fin_go_aliases = kws_all['go_aliases']
            if os.path.exists(fin_go_aliases):
                go2letter = read_d1_letter(fin_go_aliases)
                if go2letter:
                    kws_dag['go2letter'] = go2letter
        return kws_dag

    @staticmethod
    def _get_tcntobj(goids, go2obj, **kws):
        """Get a TermCounts object if the user provides an annotation file, otherwise None."""
        # kws: gaf (gene2go taxid)
        if not AnnoReaderBase.valid_formats.isdisjoint(kws):
            # Get a reduced go2obj set for TermCounts
            _gosubdag = GoSubDag(goids, go2obj, rcntobj=False)
            kws = dict(kws)
            kws['godag'] = go2obj
            return get_tcntobj(go2obj, **kws)  # TermCounts

    def get_docargs(self, args=None, prt=None):
        """Pare down docopt. Return a minimal dictionary and a set containing runtime arg values."""
        docargs = self.objdoc.get_docargs(args, prt)
        self._chk_docopts(docargs)
        return docargs

    def _chk_docopts(self, kws):
        """Check for common user command-line errors."""
        # outfile should contain .png, .png, etc.
        outfile = kws['outfile']
        if len(kws) == 2 and os.path.basename(kws['obo']) == "go-basic.obo" and \
            kws['outfile'] == self.dflt_outfile:
            self._err("NO GO IDS SPECFIED", err=False)
        if 'obo' in outfile:
            self._err("BAD outfile({O})".format(O=outfile))
        if 'gaf' in kws and 'gene2go' in kws:
            self._err("SPECIFY ANNOTAIONS FROM ONE FILE")
        if 'gene2go' in kws:
            if 'taxid' not in kws:
                self._err("SPECIFIY taxid WHEN READ NCBI'S gene2go FILE")

    def _err(self, msg, err=True):
        """Print useage and error before exiting."""
        severity = "FATAL" if err else "NOTE"
        txt = "".join([self.objdoc.doc,
                       "User's command-line:\n\n",
                       "  % go_plot.py {ARGS}\n\n".format(ARGS=" ".join(sys.argv[1:])),
                       "**{SEV}: {MSG}\n".format(SEV=severity, MSG=msg)])
        if err:
            raise RuntimeError(txt)
        sys.stdout.write(txt)
        sys.exit(0)

    def get_outfile(self, outfile, goids=None):
        """Return output file for GO Term plot."""
        # 1. Use the user-specfied output filename for the GO Term plot
        if outfile is not None and outfile != self.dflt_outfile:
            return outfile
        # 2. If only plotting 1 GO term, use GO is in plot name
        if goids is not None and len(goids) == 1:
            goid = next(iter(goids))
            goobj = self.gosubdag.go2obj[goid]
            fout = "GO_{NN}_{NM}".format(NN=goid.replace("GO:", ""), NM=goobj.name)
            return ".".join([re.sub(r"[\s#'()+,-./:<=>\[\]_}]", '_', fout), 'png'])
        # 3. Return default name
        return self.dflt_outfile

    @staticmethod
    def _get_optional_attrs(kws):
        """Given keyword args, return optional_attributes to be loaded into the GODag."""
        # Ex: def defn synonym relationship xref subset comment
        vals = OboOptionalAttrs.attributes.intersection(kws.keys())
        if 'relationships' in kws:
            vals.add('relationship')
        if 'sections' in kws:
            vals.add('relationship')
        if 'norel' in kws:
            vals.discard('relationship')
        return vals

    @staticmethod
    def _get_relationships(kws_all, relationship_in_godag):
        """Return value for GoSubDag arg, relationships."""
        if not relationship_in_godag:
            return None
        if 'relationship' in kws_all:
            return RELATIONSHIP_SET
        if 'relationships' not in kws_all:
            return None
        relationships_arg = kws_all['relationships']
        if isinstance(relationships_arg, str):
            relationships = set(kws_all['relationships'].split(','))
            chk_relationships(relationships)
            return relationships
        if relationships_arg:
            return True


# Copyright (C) 2016-2019, DV Klopfenstein, H Tang. All rights reserved.
