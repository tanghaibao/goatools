"""Product gene lists with ASCII art sections and GO IDs for each gene product."""

import sys
import collections as cx
from goatools.utils import get_b2aset
from goatools.rpt.goea_nt_xfrm import MgrNtGOEAs
from goatools.grouper.grprobj import Grouper
from goatools.grouper.sorter import Sorter
from goatools.grouper.wrxlsx import WrXlsxSortedGos

__copyright__ = "Copyright (C) 2016-2018, DV Klopfenstein, H Tang, All rights reserved."
__author__ = "DV Klopfenstein"


# pylint: disable=too-many-instance-attributes
class AArtGeneProductSetsOne(object):
    """Product gene lists with ASCII art sections and GO IDs for each gene product."""
    # nts need: nt.GO and nt.study_items

    def __init__(self, name, goea_results, obj):
        self.name = name
        self.datobj = obj  # AArtGeneProductSetsAll
        _ini = _Init(obj)
        self.go2nt = _ini.get_go2nt(goea_results)
        _grprobj = Grouper("grp", self.go2nt, obj.hdrobj, obj.grprdflt.gosubdag, go2nt=self.go2nt)
        self.sortobj = Sorter(_grprobj)
        self.sec2gos = _ini.get_sec2gos(self.sortobj)
        self.sec2chr = cx.OrderedDict([(s, obj.sec2chr[s]) for s in self.sec2gos.keys()])
        self.go2chrs = _ini.get_go2chrs(self.sec2gos, self.sec2chr)
        self.gene2gos = _ini.get_gene2gos(self.go2nt)
        self.gene2section2gos = _ini.get_gene2section2gos(self.gene2gos, self.sec2gos)
        self.gene2aart = _ini.get_gene2aart(self.gene2section2gos, self.sec2chr)

    def prt_report_grp0(self, prt=sys.stdout):
        """Print full GO/gene report without grouping."""
        summaryline = self.str_summaryline()
        kws_grp = {'use_sections':False,
                   'hdrgo_prt':False,
                   'sortby':lambda nt: [-1*nt.dcnt, nt.depth]}
        # Print grouped GO IDs
        prt.write("{SUMMARY}\n".format(SUMMARY=summaryline))
        self.prt_gos_grouped(sys.stdout, **kws_grp)
        # genes
        genes = sorted(self.gene2gos.keys())
        prt.write("\n\n{SUMMARY}\n\n".format(SUMMARY=summaryline))
        self.prt_gene_aart(genes, prt)
        # Sort genes
        prt.write("\n\n{SUMMARY}\n\n".format(SUMMARY=summaryline))
        self.prt_gene_aart_details(genes, prt)
        return (self.name, self.get_section_marks())

    def prt_report_grp1(self, prt=sys.stdout, **kws_grp):
        """Print full GO/gene report with grouping."""
        summaryline = self.str_summaryline()
        # Print grouped GO IDs
        prt.write("{SUMMARY}\n".format(SUMMARY=summaryline))
        self.prt_gos_grouped(prt, **kws_grp)
        # genes
        genes = sorted(self.gene2gos.keys())
        prt.write("\n\n{SUMMARY}\n\n".format(SUMMARY=summaryline))
        self.prt_section_key(prt)
        self.prt_gene_aart(genes, prt)
        # Sort genes
        prt.write("\n\n{SUMMARY}\n\n".format(SUMMARY=summaryline))
        self.prt_gene_aart_details(genes, prt)
        return (self.name, self.get_section_marks())

    def str_summaryline(self):
        """Print: 47 GOs, 262 genes described by 10 of 19 sections consistent_increase."""
        return "{N} GOs, {M} genes described by {X} of {Y} sections {NM}".format(
            N=len(self.go2nt), M=len(self.gene2gos),
            X=len(self.sec2chr), Y=len(self.datobj.sec2chr), NM=self.name)

    def prt_gos_grouped(self, prt, **kws_grp):
        """Print grouped GO list."""
        prtfmt = self.datobj.kws['fmtgo']
        wrobj = WrXlsxSortedGos(self.name, self.sortobj)
        # Keyword arguments: control content: hdrgo_prt section_prt top_n use_sections
        desc2nts = self.sortobj.get_desc2nts(**kws_grp)
        wrobj.prt_txt_desc2nts(prt, desc2nts, prtfmt)

    def prt_gos_flat(self, prt):
        """Print flat GO list."""
        prtfmt = self.datobj.kws['fmtgo']
        _go2nt = self.sortobj.grprobj.go2nt
        go2nt = {go:_go2nt[go] for go in self.go2nt}
        prt.write("\n{N} GO IDs:\n".format(N=len(go2nt)))
        _sortby = self._get_sortgo()
        for ntgo in sorted(go2nt.values(), key=_sortby):
            prt.write(prtfmt.format(**ntgo._asdict()))
        #print("FFFMMMTTT", prtfmt)

    def _get_sortgo(self):
        """Get function for sorting GO terms in a list of namedtuples."""
        if 'sortgo' in self.datobj.kws:
            return self.datobj.kws['sortgo']
        return self.datobj.grprdflt.gosubdag.prt_attr['sort'] + "\n"

    def prt_gene_aart(self, geneids, prt=sys.stdout):
        """For each gene, print ASCII art which represents its associated GO IDs."""
        patgene = self.datobj.kws["fmtgene"]
        itemid2name = self.datobj.kws.get("itemid2name")
        prt.write("\n{HDR}\n".format(HDR=self.str_hdr()))
        for geneid in geneids:
            symbol = "" if itemid2name is None else itemid2name.get(geneid, "")
            prt.write(patgene.format(AART=self.gene2aart[geneid], ID=geneid, NAME=symbol))

    def prt_gene_aart_details(self, geneids, prt=sys.stdout):
        """For each gene, print ASCII art which represents its associated GO IDs."""
        _go2nt = self.sortobj.grprobj.go2nt
        patgene = self.datobj.kws["fmtgene2"]
        patgo = self.datobj.kws["fmtgo2"]
        itemid2name = self.datobj.kws.get("itemid2name")
        chr2i = self.datobj.get_chr2idx()
        for geneid in geneids:
            gos_gene = self.gene2gos[geneid]
            symbol = "" if itemid2name is None else itemid2name.get(geneid, "")
            prt.write("\n")
            prt.write(patgene.format(AART=self.gene2aart[geneid], ID=geneid, NAME=symbol))
            go2nt = {go:(_go2nt[go], "".join(self.go2chrs[go])) for go in gos_gene}
            for ntgo, abc in sorted(go2nt.values(),
                                    key=lambda t: [chr2i[t[1][:1]], t[0].NS, -1*t[0].dcnt]):
                prt.write("{ABC} ".format(ABC=abc))
                prt.write(patgo.format(**ntgo._asdict()))

    def prt_section_key(self, prt=sys.stdout):
        """Print the section name and its alias."""
        for section_name, letter in self.datobj.sec2chr.items():
            mrk = '*' if section_name in self.sec2chr else ""
            prt.write("{M:1} {ABC} {SECT}\n".format(M=mrk, ABC=letter, SECT=section_name))

    def str_hdr(self):
        """Return a string representing the section headers: """
        return "".join([c for _, c in self.sec2chr.items()])

    def get_section_marks(self):
        """For each section in AArtGeneProducts, return '*' or "" ."""
        return [abc if s in self.sec2chr else "." for s, abc in self.datobj.sec2chr.items()]

    def get_gene2binvec(self):
        """Return a boolean vector for each gene representing GO section membership."""
        _sec2chr = self.sec2chr
        return {g:[s in s2gos for s in _sec2chr] for g, s2gos in self.gene2section2gos.items()}


class _Init(object):

    def __init__(self, objaartall):
        self.objaartall = objaartall  # AArtGeneProductSetsAll

    def get_go2nt(self, goea_results):
        """Return go2nt with added formatted string versions of the P-values."""
        go2obj = self.objaartall.grprdflt.gosubdag.go2obj
        # Add string version of P-values
        goea_nts = MgrNtGOEAs(goea_results).get_nts_strpval()
        return {go2obj[nt.GO].id:nt for nt in goea_nts if nt.GO in go2obj}

    @staticmethod
    def get_sec2gos(sortobj):
        """Initialize section_name2goids."""
        sec_gos = []
        for section_name, nts in sortobj.get_desc2nts_fnc(hdrgo_prt=True)['sections']:
            sec_gos.append((section_name, set(nt.GO for nt in nts)))
        return cx.OrderedDict(sec_gos)

    @staticmethod
    def get_gene2gos(go2nt):
        """Create a gene product to GO set dict."""
        gene2gos = cx.defaultdict(set)
        nt0 = next(iter(go2nt.values()))
        b_str = isinstance(nt0.study_items, str)
        # print("NNNNTTTT", nt0)
        for goid, ntgo in go2nt.items():
            study_items = ntgo.study_items.split(', ') if b_str else ntgo.study_items
            for geneid in study_items:
                gene2gos[geneid].add(goid)
        if b_str:
          b_set = set(isinstance(g.isdigit(), int) for g in nt0.study_items.split(', '))
          if b_set == set([True]):
            return {int(g):gos for g, gos in gene2gos.items()}
        return {g:gos for g, gos in gene2gos.items()}

    @staticmethod
    def get_go2chrs(sec2gos, sec2chr):
        """Dict: given a GO return a set of letters representing it's section membership(s)."""
        go2chrs = {}
        for goid, sections in get_b2aset(sec2gos).items():
            go2chrs[goid] = set(sec2chr[s] for s in sections)
        return go2chrs

    @staticmethod
    def get_gene2aart(gene2section2gos, sec2chr):
        """Return a string for each gene representing GO section membership."""
        geneid2str = {}
        for geneid, section2gos_gene in gene2section2gos.items():
            letters = [abc if s in section2gos_gene else "." for s, abc in sec2chr.items()]
            geneid2str[geneid] = "".join(letters)
        return geneid2str

    @staticmethod
    def get_gene2section2gos(gene2gos, sec2gos):
        """Get a list of section aliases for each gene product ID."""
        gene2section2gos = {}
        for geneid, gos_gene in gene2gos.items():
            section2gos = {}
            for section_name, gos_sec in sec2gos.items():
                gos_secgene = gos_gene.intersection(gos_sec)
                if gos_secgene:
                    section2gos[section_name] = gos_secgene
            gene2section2gos[geneid] = section2gos
        return gene2section2gos


# Copyright (C) 2016-2018, DV Klopfenstein, H Tang, All rights reserved.
