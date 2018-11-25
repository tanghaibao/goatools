"""Prints a Python sections file."""

import sys
from goatools.wr_tbl import prt_txt
from goatools.grouper.tasks import SummarySec2dHdrGos

__copyright__ = "Copyright (C) 2016-2018, DV Klopfenstein, H Tang, All rights reserved."
__author__ = "DV Klopfenstein"



class WrSectionsBase(object):
    """Tasks for writing a sections file."""

    def __init__(self, grprobj, ver_list=None):
        self.ver_list = ver_list
        self.grprobj = grprobj
        self.gosubdag = grprobj.gosubdag
        self.fncsortnt = self._init_fncsortnt(self.gosubdag.prt_attr['flds'])
        self.prtfmt = self._init_prtfmt("fmta")

    def prt_ver(self, prt):
        """Print version of GO-DAG for the GO and for GO slims."""
        if self.ver_list is not None:
            prt.write("# Versions:\n#    {VER}\n\n".format(VER="\n#    ".join(self.ver_list)))

    def get_sections_2dnt(self, sec2d_go):
        """Return a sections list containing sorted lists of namedtuples."""
        return [(nm, self.get_ntgos_sorted(gos)) for nm, gos in sec2d_go]

    def get_ntgos_sorted(self, hdrgos):
        """Return sorted Grouper namedtuples if there are user GO IDs underneath."""
        go2nt = self.grprobj.go2nt
        return sorted([go2nt[go] for go in hdrgos if go in go2nt], key=self.fncsortnt)

    def prt_ntgos(self, prt, ntgos):
        """Print the Grouper namedtuples."""
        for ntgo in ntgos:
            key2val = ntgo._asdict()
            prt.write("{GO_LINE}\n".format(GO_LINE=self.prtfmt.format(**key2val)))

    def _init_prtfmt(self, key="fmta"):
        """Return print format for Grouper, which includes hdr1usr01 and num_usrgos."""
        prtfmt = self.gosubdag.prt_attr[key]
        return prtfmt.replace("{NS}", "{NS} {hdr1usr01:2} {num_usrgos:>4} uGOs")

    def get_summary_data(self, sec2d_nt):
        """Get placed/unplaced GO IDs and sections."""
        grouped = set()
        ungrouped = set()
        sections = set()
        secdflt = self.grprobj.hdrobj.secdflt
        for section_name, nts in sec2d_nt:
            if section_name != secdflt:
                if nts:
                    grouped.update(set(nt.GO for nt in nts))
                    sections.add(section_name)
                else:
                    ungrouped.update(set(nt.GO for nt in nts))
        return {'grouped':grouped, 'ungrouped':ungrouped, 'sections':sections}

    def get_summary_str(self, sec2d_nt):
        """Get string describing counts of placed/unplaced GO IDs and count of sections."""
        data = self.get_summary_data(sec2d_nt)
        return "{M} GO IDs placed into {N} sections; {U} unplaced GO IDs".format(
            N=len(data['sections']), M=len(data['grouped']), U=len(data['ungrouped']))

    @staticmethod
    def _init_fncsortnt(flds):
        """Return a sort function for sorting header GO IDs found in sections."""
        if 'tinfo' in flds:
            if 'D1' in flds:
                return lambda ntgo: [ntgo.NS, -1*ntgo.tinfo, ntgo.depth, ntgo.D1, ntgo.alt]
            else:
                return lambda ntgo: [ntgo.NS, -1*ntgo.tinfo, ntgo.depth, ntgo.alt]
        if 'dcnt' in flds:
            if 'D1' in flds:
                return lambda ntgo: [ntgo.NS, -1*ntgo.dcnt, ntgo.depth, ntgo.D1, ntgo.alt]
            else:
                return lambda ntgo: [ntgo.NS, -1*ntgo.dcnt, ntgo.depth, ntgo.alt]
        else:
            return lambda ntgo: [ntgo.NS, -1*ntgo.depth, ntgo.alt]


class WrSectionsPy(WrSectionsBase):
    """Holds formatting information for printing sections into a Python file."""

    def __init__(self, grprobj, ver_list=None):
        super(WrSectionsPy, self).__init__(grprobj, ver_list)
        self.prtfmt = self.prtfmt.replace('{GO}', '        "{GO}", ')

    def wr_py_sections_new(self, fout_py, doc=None):
        """Write the first sections file."""
        sections = self.grprobj.get_sections_2d()
        return self.wr_py_sections(fout_py, sections, doc)

    def wr_py_sections(self, fout_py, sections, doc=None):
        """Write sections 2-D list into a Python format list."""
        if sections is None:
            sections = self.grprobj.get_sections_2d()
        sec2d_nt = self.get_sections_2dnt(sections)  # lists of GO Grouper namedtuples
        with open(fout_py, 'w') as prt:
            self._prt_py_sections(sec2d_nt, prt, doc)
            dat = SummarySec2dHdrGos().summarize_sec2hdrgos(sections)
            sys.stdout.write(self.grprobj.fmtsum.format(
                GO_DESC='hdr', SECs=len(dat['S']), GOs=len(dat['G']),
                UNGRP=len(dat['U']), undesc="unused",
                ACTION="WROTE:", FILE=fout_py))

    def _prt_py_sections(self, sec2d_nt, prt=sys.stdout, doc=None):
        """Print sections 2-D list into a Python format list."""
        if doc is None:
            doc = 'Sections variable'
        prt.write('"""{DOC}"""\n\n'.format(DOC=doc))
        self.prt_ver(prt)
        prt.write("# pylint: disable=line-too-long\n")
        strcnt = self.get_summary_str(sec2d_nt)
        prt.write("SECTIONS = [ # {CNTS}\n".format(CNTS=strcnt))
        prt.write('    # ("New Section", [\n')
        prt.write('    # ]),\n')
        for section_name, nthdrgos in sec2d_nt:
            self._prt_py_section(prt, section_name, nthdrgos)
        prt.write("]\n")

    def _prt_py_section(self, prt, section_name, ntgos):
        """Print one section and its GO headers."""
        prt.write('    ("{SEC}", [ # {N} GO-headers\n'.format(SEC=section_name, N=len(ntgos)))
        self.prt_ntgos(prt, ntgos)
        prt.write("    ]),\n")


class WrSectionsTxt(WrSectionsBase):
    """Manages GO group headers and optionally sections containing GO group headers."""

    def __init__(self, grprobj, ver_list=None):
        super(WrSectionsTxt, self).__init__(grprobj, ver_list)

    @staticmethod
    def prt_sections(prt, sections, prtfmt, secspc=False):
        """Print GO namedtuples in their sections."""
        num_goids = 0
        for section, nts_flat in sections:
            # Add an empty line between sections, if desired
            if secspc:
                prt.write("\n")
            prt.write("{SECTION}\n".format(SECTION=section))
            num_nts = len(nts_flat)
            num_goids += num_nts
            prt_txt(prt, nts_flat, prtfmt=prtfmt)

    def prt_info(self, prt, sections=None):
        """Print GO namedtuples in their sections."""
        if sections is None:
            sections = self.grprobj.get_sections_2d()
        num_goids = 0
        for section_name, nts_flat in sections:
            num_nts = len(nts_flat)
            num_goids += num_nts
            prt.write("{N:3} GO IDs in section({SEC})\n".format(N=num_nts, SEC=section_name))
        prt.write("{N:3} GO IDs\n".format(N=num_goids))

    def prt_goid_cnt(self, prt=sys.stdout):
        """Get number of hdrgos and usrgos in each section."""
        for section_name, hdrgos_sec in self.grprobj.get_sections_2d():
            prt.write("{NAME} {Us:5,} {Hs:5,} {SEC}\n".format(
                NAME=self.grprobj.grpname,
                Us=len(self.grprobj.get_usrgos_g_hdrgos(hdrgos_sec)),
                Hs=len(hdrgos_sec),
                SEC=section_name))

    def wr_txt_grouping_gos(self):
        """Write one file per GO group."""
        prt_goids = self.grprobj.gosubdag.prt_goids
        for hdrgo, usrgos in self.grprobj.hdrgo2usrgos.items():
            keygos = usrgos.union([hdrgo])
            fout_txt = "{BASE}.txt".format(BASE=self.grprobj.get_fout_base(hdrgo))
            with open(fout_txt, 'w') as prt:
                prt_goids(keygos, prt=prt)
                sys.stdout.write("  {N:5,} GO IDs WROTE: {TXT}\n".format(
                    N=len(keygos), TXT=fout_txt))

    def wr_txt_section_hdrgos(self, fout_txt, sortby=None, prt_section=True):
        """Write high GO IDs that are actually used to group current set of GO IDs."""
        sec2d_go = self.grprobj.get_sections_2d()    # lists of GO IDs
        sec2d_nt = self.get_sections_2dnt(sec2d_go)  # lists of GO Grouper namedtuples
        if sortby is None:
            sortby = self.fncsortnt
        with open(fout_txt, 'w') as prt:
            self.prt_ver(prt)
            prt.write("# GROUP NAME: {NAME}\n".format(NAME=self.grprobj.grpname))
            for section_name, nthdrgos_actual in sec2d_nt:
                if prt_section:
                    prt.write("# SECTION: {SECTION}\n".format(SECTION=section_name))
                self.prt_ntgos(prt, nthdrgos_actual)
                if prt_section:
                    prt.write("\n")
            dat = SummarySec2dHdrGos().summarize_sec2hdrgos(sec2d_go)
            sys.stdout.write(self.grprobj.fmtsum.format(
                GO_DESC='hdr', SECs=len(dat['S']), GOs=len(dat['G']),
                UNGRP=len(dat['U']), undesc="unused",
                ACTION="WROTE:", FILE=fout_txt))
        return sec2d_nt


# Copyright (C) 2016-2018, DV Klopfenstein, H Tang, All rights reserved.
