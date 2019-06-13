"""Print a GO term's lower-level hierarchy."""

__copyright__ = "Copyright (C) 2016-2019, DV Klopfenstein, H Tang. All rights reserved."
__author__ = "DV Klopfenstein"

import sys
import collections as cx
from goatools.godag.consts import NAMESPACE2GO
from goatools.gosubdag.go_paths import GoPaths
from goatools.rpt.write_hierarchy_base import WrHierPrt


class WrHierGO(object):
    """Write hierarchy object."""

    kws_dct = set(['max_indent'])
    kws_set = set(['no_indent', 'concise'])

    def __init__(self, gosubdag, **kws):
        self.gosubdag = gosubdag  # GoSubDag arg, children=True, must be used
        self.usrdct = {k:v for k, v in kws.items() if k in kws}
        self.usrset = set([k for k, v in kws.items() if k in kws and v])
        # ' {NS} {dcnt:6,} L{level:02} D{depth:02} {D1:5} {GO_name}'

    def prt_hier_all(self, prt=sys.stdout):
        """Write hierarchy for all GO Terms in obo file."""
        # Print: [biological_process, molecular_function, and cellular_component]
        items_list = set()
        for goid in ['GO:0008150', 'GO:0003674', 'GO:0005575']:
            items_list.update(self.prt_hier_down(goid, prt))
        return items_list

    def prt_hier_down(self, goid, prt=sys.stdout):
        """Write hierarchy for all GO IDs below GO ID in arg, goid."""
        wrhiercfg = self._get_wrhiercfg()
        obj = WrHierPrt(self.gosubdag.go2obj, self.gosubdag.go2nt, wrhiercfg, prt)
        obj.prt_hier_rec(goid)
        return obj.items_list

    def prt_hier_up(self, goids, prt=sys.stdout):
        """Write hierarchy for all GO IDs below GO ID in arg, goid."""
        go2goterm_all = {go:self.gosubdag.go2obj[go] for go in goids}
        objp = GoPaths()
        items_list = []
        for namespace, go2term_ns in self._get_namespace2go2term(go2goterm_all).items():
            goids_all = set()  # GO IDs from user-specfied GO to root
            for goid_usr, goterm in go2term_ns.items():
                goids_all.add(goid_usr)
                paths = objp.get_paths_from_to(goterm, goid_end=None, dn0_up1=True)
                goids_all.update(set(o.id for p in paths for o in p))
            # Only include GO IDs from user-specified GO to the root
            if 'include_only' not in self.usrdct:
                self.usrdct['include_only'] = set()
            self.usrdct['include_only'].update(goids_all)
            # Mark the user-specfied GO term
            if 'item_marks' not in self.usrdct:
                self.usrdct['item_marks'] = {}
            for goid_usr in go2term_ns.keys():
                if goid_usr not in self.usrdct['item_marks']:
                    self.usrdct['item_marks'][goid_usr] = '*'
            # Write the hierarchy
            wrhiercfg = self._get_wrhiercfg()
            obj = WrHierPrt(self.gosubdag.go2obj, self.gosubdag.go2nt, wrhiercfg, prt)
            go_root = self._get_goroot(goids_all, namespace)
            obj.prt_hier_rec(go_root)
            items_list.extend(obj.items_list)
        return items_list

    @staticmethod
    def _get_namespace2go2term(go2terms):
        """Group GO IDs by namespace."""
        namespace2go2term = cx.defaultdict(dict)
        for goid, goterm in go2terms.items():
            namespace2go2term[goterm.namespace][goid] = goterm
        return namespace2go2term

    def _get_wrhiercfg(self):
        """Initialize print format."""
        prtfmt = self.gosubdag.prt_attr['fmt']
        prtfmt = prtfmt.replace('{GO} # ', '')
        prtfmt = prtfmt.replace('{D1:5} ', '')
        return {'name2prtfmt':{'ITEM':prtfmt, 'ID':'{GO}{alt:1}'},
                'max_indent': self.usrdct.get('max_indent'),
                'include_only': self.usrdct.get('include_only'),
                'item_marks': self.usrdct.get('item_marks', {}),
                'concise_prt': 'concise' in self.usrset,
                'indent': 'no_indent' not in self.usrset,
                'dash_len': self.usrdct.get('dash_len', 6),
                'sortby': self.usrdct.get('sortby')
               }

    def _get_goroot(self, goids_all, namespace):
        """Get the top GO for the set of goids_all."""
        root_goid = NAMESPACE2GO[namespace]
        if root_goid in goids_all:
            return root_goid
        root_goids = set()
        for goid in goids_all:
            goterm = self.gosubdag.go2obj[goid]
            if goterm.depth == 0:
                root_goids.add(goterm.id)
        if len(root_goids) == 1:
            return next(iter(root_goids))
        raise RuntimeError("UNEXPECTED NUMBER OF ROOTS: {R}".format(R=root_goids))

#### Examples:
####
#### Print the GO IDs associated with human genes
#### >>> python scripts/wr_hier.py BP --gene2go=gene2go --taxid=9606 --concise -o BP_9606.txt
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
#### >>> python {SCR} --max_indent=2 --dash_len=2
####
####
#### Print a conciseened version of the hierarchy for BP, MF, and CC.
#### This will only print a path to a leaf GO Term once.
#### If the path appears a second time, the term is printed again, but its path is not.
#### The presence of a compressed (unprinted) paths is marked by using '=" instead of '-'.
####
####     $ wc -l hier_BP_MF_CC*.rpt
####
####           789583 hier_BP_MF_CC.rpt
####            70152 hier_BP_MF_CC_concise.rpt
####
#### >>> python {SCR} --o=hier_BP_MF_CC_concise.rpt --concise
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

# Copyright (C) 2016-2019, DV Klopfenstein, H Tang. All rights reserved.
