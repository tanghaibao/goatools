"""Create a pydot Node for a GO Term."""

__copyright__ = "Copyright (C) 2016-2020, DV Klopfenstein, H Tang, All rights reserved."
__author__ = "DV Klopfenstein"

import pydot
from goatools.gosubdag.plot.go_name_shorten import ShortenText
from goatools.gosubdag.utils import extract_kwargs

class GoNodeOpts:
    """Processes GO Node plot args."""

    exp_keys = set(['goobj2fncname', 'go2txt', 'objgoea', 'prt_flds'])

    exp_elems = set([
        'c2ps',         # Count of an object's Parent
        'prt_pcnt',     # Always print parent count: pN
        'parentcnt',    # Print parent count only if not all parents are shown
        'childcnt',     # Always print child count: cN
        'mark_alt_id',  # Put an 'a' after GO:NNNNNNN if it is an alternate GO ID
        'shorten',      # Shorten GO description
        'no_name',      # Do not print GO description
    ])

    def __init__(self, gosubdag, **kws):
        self.gosubdag = gosubdag
        # kws = {'set':set(...), 'dict':{...}}
        self.kws = extract_kwargs(kws, self.exp_keys, self.exp_elems)

    def get_kws(self):
        """Only load keywords if they are specified by the user."""
        ret = self.kws['dict'].copy()
        act_set = self.kws['set']
        if 'shorten' in act_set and 'goobj2fncname' not in ret:
            ret['goobj2fncname'] = ShortenText().get_short_plot_name
        if 'dict' in self.kws and 'go2txt' in self.kws['dict']:
            self._init_go2txt_altgos(self.kws['dict']['go2txt'])
        return ret

    def get_present(self):
        """Only store keywords if they are specified by the user."""
        # The presence of c2ps marks that the user specified parentcnt=True
        return self.kws['set'].difference(['parentcnt'])

    def _init_go2txt_altgos(self, go2txt):
        """If user provided GO.alt_id, add the corressponding main GO ID, if needed"""
        _go2obj = self.gosubdag.go2obj
        for goid_user, txt in go2txt.items():
            if goid_user in _go2obj:
                goid_main = _go2obj[goid_user].item_id
                if goid_user != goid_main and goid_main not in go2txt:
                    go2txt[goid_main] = txt


class GoNode:
    """Creates pydot Node containing a GO term."""

    exclude = {'tfreq',}

    def __init__(self, gosubdag, objcolor, optobj):
        self.gosubdag = gosubdag     # GoSubDag
        self.objcolor = objcolor     # Go2Color   -> color options
        self.kws = optobj.get_kws()  # GoNodeOpts -> text  options
        self.present = optobj.get_present()
        self.go2color = objcolor.go2color

    def get_node(self, goid, goobj):
        """Return pydot node."""
        # pydot.Node.objdict holds this information. pydot.Node.objdict['name']
        return pydot.Node(
            self.get_node_text(goid, goobj),
            shape="box",
            style="rounded, filled",
            fillcolor=self.go2color.get(goid, "white"),
            color=self.objcolor.get_bordercolor(goid))

    def str_fmthdr(self, goid, goobj):
        """Return hdr line seen inside a GO Term box."""
        # Shorten: Ex: GO:0007608 -> G0007608
        go_txt = goid.replace("GO:", "G")
        if 'mark_alt_id' in self.present and goid != goobj.id:
            go_txt += 'a'
        return go_txt

    # ----------------------------------------------------------------------------------
    # Methods for text printed inside GO terms
    def get_node_text(self, goid, goobj):
        """Return a string to be printed in a GO term box."""
        txt = []
        # Header line: "GO:0036464 L04 D06"
        hdr = self.get_hdr(goid, goobj)
        if hdr != '':
            txt.append(hdr)
        # GO name line: "cytoplamic ribonucleoprotein"
        if 'no_name' not in self.present:
            txt.append(self._get_go_name(goobj))
        # study info line: "24 genes"
        if 'objgoea' in self.kws:
            study_txt = self.kws['objgoea'].get_study_txt(goid)
            if study_txt is not None:
                txt.append(study_txt)
        # Add user-specified text, if needed
        if 'go2txt' in self.kws and goid in self.kws['go2txt']:
            txt.append(self.kws['go2txt'][goid])
        return "\n".join(txt)

    def _get_go_name(self, goobj):
        """Return GO name/description, as is or edited by a user function."""
        if 'goobj2fncname' not in self.kws:
            return goobj.name.replace(",", "\n")
        # Return GO Term name edited by user-provided function
        return self.kws['goobj2fncname'](goobj)

    def get_hdr(self, goid, goobj):
        """Header for GO Term box. Ex: 'G0001719 L6 D9 d3.'"""
        hdr = []
        ntgo = self.gosubdag.go2nt.get(goid)
        prt_flds = self._get_prtflds()
        # Add letter to depth-01 GO Node.
        if 'D1' in prt_flds and goobj.depth == 1:
            hdr.append("{ABC} ".format(ABC=ntgo.D1))
        if 'GO' in prt_flds:
            hdr.append(self.str_fmthdr(goid, goobj))
        if 'level' in prt_flds:
            hdr.append("L{level}".format(level=goobj.level))
        if 'depth' in prt_flds:
            hdr.append("D{depth}".format(depth=goobj.depth))
        if 'reldepth' in prt_flds:
            hdr.append("R{reldepth}".format(reldepth=goobj.reldepth))
        # Print count of parents for this GO term
        if 'c2ps' in self.kws:
            self._add_parent_cnt(hdr, goobj, self.kws['c2ps'])
        # Print count of children for this GO term
        childcnt_str = self._get_hdr_childcnt(goobj, ntgo)
        if childcnt_str:
            hdr.append(childcnt_str)
        # Print count of all descendants down to the leaf-level for this GO term
        if 'dcnt' in prt_flds:
            hdr.append("d{N}".format(N=ntgo.dcnt))
        if 'tinfo' in prt_flds:
            hdr.append("i{I:4.02f}".format(I=ntgo.tinfo))
        if 'tfreq' in prt_flds:
            hdr.append("f{I:4.03f}".format(I=ntgo.tfreq))
        if 'REL' in prt_flds:
            hdr.append("{R}".format(R=ntgo.REL_short))
        return " ".join(hdr) if hdr else ''

    def _get_prtflds(self):
        """Get print fields for GO header."""
        # User-specified print fields
        ntflds = self.gosubdag.prt_attr['flds']
        prt_flds = self.kws.get('prt_flds')
        if prt_flds:
            return prt_flds.intersection(ntflds)
        # Default print fields
        exclude = set(self.exclude)
        if self.gosubdag.relationships:
            exclude.add('level')
        return set(f for f in ntflds if f not in exclude)

    def _get_hdr_childcnt(self, goobj, ntgo):
        """Get string representing count of children for this GO term."""
        if 'childcnt' in self.present:
            return "c{N}".format(N=len(goobj.children))
        if self.gosubdag.relationships and not goobj.children and ntgo.dcnt != 0:
            return "c0"
        return None

    def _add_parent_cnt(self, hdr, goobj, c2ps):
        """Add the parent count to the GO term box for if not all parents are plotted."""
        if goobj.id in c2ps:
            parents = c2ps[goobj.id]
            if 'prt_pcnt' in self.present or parents and len(goobj.parents) != len(parents):
                assert len(goobj.parents) == len(set(goobj.parents))
                hdr.append("p{N}".format(N=len(set(goobj.parents))))


# Copyright (C) 2016-2020, DV Klopfenstein, H Tang, All rights reserved.
