"""Print a GO term's lower-level hierarchy."""

__copyright__ = "Copyright (C) 2016-2018, DV Klopfenstein, H Tang. All rights reserved."
__author__ = "DV Klopfenstein"

import sys
import collections as cx
from goatools.godag.consts import Consts
from goatools.gosubdag.go_paths import GoPaths


class WrHierGO(object):
    """Write hierarchy object."""

    kws_dct = set(['max_indent', 'num_child'])
    kws_set = set(['no_indent', 'short'])
    consts = Consts()

    def __init__(self, gosubdag, **kws):
        self.gosubdag = gosubdag  # GoSubDag arg, children=True, must be used
        self.usrdct = {k:v for k, v in kws.items() if k in kws}
        self.usrset = set([k for k, v in kws.items() if k in kws and v])
        # ' {NS} {dcnt:6,} L{level:02} D{depth:02} {D1:5} {GO_name}'

    def prt_hier_all(self, prt=sys.stdout):
        """Write hierarchy for all GO Terms in obo file."""
        # Print: [biological_process, molecular_function, and cellular_component]
        for goid in ['GO:0008150', 'GO:0003674', 'GO:0005575']:
            self.prt_hier_down(goid, prt)

    def prt_hier_down(self, goid, prt=sys.stdout):
        """Write hierarchy for all GO IDs below GO ID in arg, goid."""
        obj = _WrHierPrt(self, prt)
        obj.prt_hier_rec(goid)

    def prt_hier_up(self, goids, prt=sys.stdout):
        """Write hierarchy for all GO IDs below GO ID in arg, goid."""
        go2goterm_all = {go:self.gosubdag.go2obj[go] for go in goids}
        objp = GoPaths()
        for namespace, go2term_ns in self._get_namespace2go2term(go2goterm_all).items():
            go_root = self.consts.NAMESPACE2GO[namespace]
            goids_all = set()
            for goid, goterm in go2term_ns.items():
                goids_all.add(goid)
                paths = objp.get_paths_from_to(goterm, goid_end=None, dn0_up1=True)
                goids_all.update(set(o.id for p in paths for o in p))
            obj = _WrHierPrt(self, prt, goids_all, set(go2term_ns.keys()))
            obj.prt_hier_rec(go_root)

    @staticmethod
    def _get_namespace2go2term(go2terms):
        """Group GO IDs by namespace."""
        namespace2go2term = cx.defaultdict(dict)
        for goid, goterm in go2terms.items():
            namespace2go2term[goterm.namespace][goid] = goterm
        return namespace2go2term


class _WrHierPrt(object):
    """Print GO hierarchy."""

    def __init__(self, obj, prt=sys.stdout, include_only=None, go_marks=None):
        self.gosubdag = obj.gosubdag
        self.max_indent = obj.usrdct.get('max_indent')
        # self.num_child = num_child
        self.include_only = include_only  # if include_only else set()
        self.go_marks = go_marks if go_marks else set()
        self.short_prt = 'short' in obj.usrset
        # vars
        self.prt = prt
        self.gos_printed = set()
        self.prtfmt = self.gosubdag.prt_attr['fmta'].replace('{GO}{alt:1} # ', '')

    def prt_hier_rec(self, goid, depth=1):
        """Write hierarchy for a GO Term record and all GO IDs down to the leaf level."""
        ntgo = self.gosubdag.go2nt[goid]
        ntobj = self.gosubdag.go2obj[goid]
        # Shortens hierarchy report by only printing the hierarchy
        # for the sub-set of user-specified GO terms which are connected.
        if self.include_only is not None and goid not in self.include_only:
            return
        nrp = self.short_prt and goid in self.gos_printed
        if self.go_marks:
            self.prt.write('{} '.format('>' if goid in self.go_marks else ' '))

        # '-' is default character indicating hierarchy level
        # '=' is used to indicate a hierarchical path printed in detail previously.
        dashgo = self._str_dashgoid(ntgo, depth, not nrp or not ntobj.children)
        self.prt.write('{DASHGO:{N}}'.format(DASHGO=dashgo, N=17))

        # if num_child is not None:
        #     self.prt.write('{N:>5} '.format(N=len(self.get_all_children())))
        self.prt.write("{GO_INFO}\n".format(GO_INFO=self.prtfmt.format(**ntgo._asdict())))
        self.gos_printed.add(goid)
        # Do not print hierarchy below this turn if it has already been printed
        if nrp:
            return
        depth += 1
        if self.max_indent is not None and depth > self.max_indent:
            return
        for child in ntobj.children:
            self.prt_hier_rec(child.id, depth)

    @staticmethod
    def _str_dash(depth, single_or_double):
        """Return a string containing dashes (optional) and GO ID."""
        # '-' is default character indicating hierarchy level
        # '=' is used to indicate a hierarchical path printed in detail previously.
        letter = '-' if single_or_double else '='
        return ''.join([letter]*depth)

    def _str_dashgoid(self, ntgo, depth, single_or_double):
        """Return a string containing dashes (optional) and GO ID."""
        dashes = self._str_dash(depth, single_or_double)
        return "{DASHES} {GO}{alt:1}".format(DASHES=dashes, GO=ntgo.GO, alt=ntgo.alt)


#### Examples:
####
#### Print the hierarchy below Term, GO:0030663
#### >>> python {SCR} GO:0030663
####
#### - GO:0030663	level-05	depth-07	COPI-coated vesicle membrane [cellular_component]
#### -- GO:0012508	level-05	depth-08	Golgi to ER transport vesicle membrane [cellular_component]
#### -- GO:0012509	level-05	depth-08	inter-Golgi transport vesicle membrane [cellular_component]
####
####
#### Write the hierarchy below Term, GO:0030663 into a file
#### >>> python {SCR} GO:0030663 --o=hier_GO_0030663.rpt
####
####   WROTE: hier_GO_0030663.rpt
####
#### Print the hierarchy for biological process, molecular_function, and cellular_component:
#### >>> python {SCR} --o=hier_BP_MF_CC.rpt
####
#### Print hierarchy for BP, MF, CC only printing the first 2 levels.
#### >>> python {SCR} --max_indent=2
#### >>> python {SCR} --max_indent=2 --dash_len=2 --num_child
####
####
#### Print a shortened version of the hierarchy for BP, MF, and CC.
#### This will only print a path to a leaf GO Term once.
#### If the path appears a second time, the term is printed again, but its path is not.
#### The presence of a compressed (unprinted) paths is marked by using '=" instead of '-'.
####
####     $ wc -l hier_BP_MF_CC*.rpt
####
####           789583 hier_BP_MF_CC.rpt
####            70152 hier_BP_MF_CC_short.rpt
####
#### >>> python {SCR} --o=hier_BP_MF_CC_short.rpt --short
####
#### Print hierarchy
#### -  26894 GO:0008150	level-00	depth-00	biological_process [biological_process]
#### --    30 GO:0001906	level-01	depth-01	cell killing [biological_process]
#### --   555 GO:0002376	level-01	depth-01	immune system process [biological_process]
#### -- 11208 GO:0065007	level-01	depth-01	biological regulation [biological_process]
####
#### >>> python {SCR}
####
#### This program prints the hierarchy for all GO terms, if no argument is provided.
#### If a GO term is provided as an argument, then the hierarchy of all children
#### for that term is printed.
####
#### """.format(SCR='write_hierarchy')


# Copyright (C) 2016-2018, DV Klopfenstein, H Tang. All rights reserved.
