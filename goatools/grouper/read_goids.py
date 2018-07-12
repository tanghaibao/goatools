"""Functions to read text or tsv files containing GO IDs and sections of GO IDs."""

from __future__ import print_function

import os
import sys
import re
import pkgutil
import importlib
from goatools.gosubdag.go_tasks import chk_goids
from goatools.grouper.hdrgos import HdrgosSections
from goatools.grouper.grprobj import Grouper
from goatools.grouper.tasks import SummarySec2dHdrGos

__copyright__ = "Copyright (C) 2016-2018, DV Klopfenstein, All rights reserved."
__author__ = "DV Klopfenstein"


def read_sections(sections_file, exclude_ungrouped=False, prt=sys.stdout):
    """Get sections and GO grouping hdrgos from file, if sections exist."""
    if sections_file is None:
        return None
    assert isinstance(sections_file, str), "BAD SECTIONS FILENAME({S})".format(
        S=sections_file)
    if os.path.exists(sections_file):
        return ReadGoids().read_sections(sections_file, False, exclude_ungrouped)
    # Is 'sections_file' a module string?
    if '/' not in sections_file and r'\\' not in sections_file and \
            pkgutil.find_loader(sections_file) is not None:
        mod = importlib.import_module(sections_file)
        var = getattr(mod, 'SECTIONS', None)
        if var is not None:
            dat = SummarySec2dHdrGos().summarize_sec2hdrgos(var)
            print(Grouper.fmtsum.format(
                GO_DESC='hdr', SECs=len(dat['S']), GOs=len(dat['G']),
                UNGRP="N/A", undesc="unused", ACTION="IMPORTED: ", FILE=sections_file))
            return var
        raise RuntimeError("NO 'SECTIONS' VARIABLE FOUND IN MODULE({M})".format(M=sections_file))
    if prt:
        prt.write("CANNOT READ: {SEC}\n".format(SEC=sections_file))


def read_goids(fin_txt, get_goids_only=False, exclude_ungrouped=False, prt=sys.stdout):
    """Get user list of GO IDs either from a list or from GO IDs on the command-line"""
    return ReadGoids().read_txt(fin_txt, get_goids_only, exclude_ungrouped, prt)


class ReadGoids(object):
    """Get user list of GO IDs either from a list or from GO IDs on the command-line"""

    srch_section = re.compile(r'^#?\s*SECTION:\s*(\S.*\S)\s*$', flags=re.IGNORECASE)
    # For reading SECTIONS from a Python file as text (without importing)
    #   ("cell death", [ # 6 GO-headers
    #   ("viral/bacteria", [ # 4 GO-headers
    srch_py_section = re.compile(r'^\s*(\(|\[)\s*(\'|")(.*)\s*(\'|")\s*,\s*\[')
    #       "GO:0002376", # BP    564 L01 D01 M     immune system process
    #       "GO:0002682", # BP  1,183 L02 D02 AB    regulation of immune system process
    srch_py_goids = re.compile(r'^\s*(\'|")(GO:\d{7})(\'|")\s*,?')

    def __init__(self):
        self.goids_fin = []
        self.sections_seen = []
        self.section2goids = {}

    def read_txt(self, fin_txt, get_goids_only, exclude_ungrouped, prt=sys.stdout):
        """Get user list of GO IDs either from a list or from GO IDs on the command-line"""
        goids_fin = self._read_txt(fin_txt, get_goids_only, exclude_ungrouped)
        sections = self._read_finish(goids_fin, prt)
        # Print summary of GO IDs read
        if prt is not None:
            self._prt_read_msg(prt, fin_txt, exclude_ungrouped)
        return sections

    def read_py(self, fin_txt, get_goids_only, exclude_ungrouped, prt=sys.stdout):
        """Read GO IDs or sections data from a Python file."""
        goids_fin = self._read_py(fin_txt, get_goids_only, exclude_ungrouped)
        sections = self._read_finish(goids_fin, prt)
        # Print summary of GO IDs read
        if prt is not None:
            self._prt_read_msg(prt, fin_txt, exclude_ungrouped)
        return sections

    def read_sections(self, sections_file, get_goids_only, exclude_ungrouped):
        """Read sections variable from a text file of from a Python file."""
        ext = os.path.splitext(sections_file)[1]
        file_contents = None
        if ext == '.py':
            file_contents = self.read_py(sections_file, get_goids_only, exclude_ungrouped)
        else:
            file_contents = self.read_txt(sections_file, get_goids_only, exclude_ungrouped)
        if file_contents:
            return file_contents.get('sections', None)

    def _read_finish(self, goids_fin, prt):
        """Get one of: {'goids':...} or {'sections':...} from reading a file."""
        # Report unused sections, if any
        if len(self.section2goids) != len(self.sections_seen):
            self._rpt_unused_sections(prt)
        # If there are no sections, then goids_fin holds all GO IDs in file
        if not self.sections_seen:
            self.goids_fin = goids_fin

        if goids_fin:
            return self.internal_get_goids_or_sections()  # {'goids':...} or {'sections':...}
        else:
            sys.stdout.write(
                "\n**WARNING: GO IDs MUST BE THE FIRST 10 CHARACTERS OF EACH LINE\n\n")

    def _read_txt(self, fin_txt, get_goids_only, exclude_ungrouped):
        """Read GO file. Store results in: section2goids sections_seen. Return goids_fin."""
        goids_sec = []
        with open(fin_txt) as istrm:
            # Lines starting with a GO ID will have that GO ID read and stored.
            #   * Lines that do not start with a GO ID will be ignored.
            #   * Text after the 10 characters in a GO ID will be ignored.
            section_name = None
            for line in istrm:
                if line[:3] == "GO:":
                    goids_sec.append(line[:10])
                elif not get_goids_only and ":" in line:
                    mtch = self.srch_section.match(line)
                    if mtch:
                        secstr = mtch.group(1)
                        if section_name is not None and goids_sec:
                            self.section2goids[section_name] = goids_sec
                        if not exclude_ungrouped or secstr != HdrgosSections.secdflt:
                            section_name = secstr
                            self.sections_seen.append(section_name)
                        else:
                            section_name = None
                        goids_sec = []
            if section_name is not None and goids_sec:
                self.section2goids[section_name] = goids_sec
        return goids_sec

    def _read_py(self, fin_py, get_goids_only, exclude_ungrouped):
        """Read Python sections file. Store: section2goids sections_seen. Return goids_fin."""
        goids_sec = []
        with open(fin_py) as istrm:
            section_name = None
            for line in istrm:
                mgo = self.srch_py_goids.search(line)     # Matches GO IDs in sections
                if mgo:
                    goids_sec.append(mgo.group(2))
                elif not get_goids_only and "[" in line:
                    msec = self.srch_py_section.search(line)  # Matches sections
                    if msec:
                        secstr = msec.group(3)
                        if section_name is not None and goids_sec:
                            self.section2goids[section_name] = goids_sec
                        if not exclude_ungrouped or secstr != HdrgosSections.secdflt:
                            section_name = secstr
                            self.sections_seen.append(section_name)
                        else:
                            section_name = None
                        goids_sec = []
            if section_name is not None and goids_sec:
                self.section2goids[section_name] = goids_sec
        return goids_sec

    def _rpt_unused_sections(self, prt):
        """Report unused sections."""
        sections_unused = set(self.sections_seen).difference(self.section2goids.keys())
        for sec in sections_unused:
            prt.write("  UNUSED SECTION: {SEC}\n".format(SEC=sec))

    def internal_get_goids_or_sections(self):
        """Return GO IDs, Sections/GOs, or None."""
        if self.goids_fin:
            chk_goids(self.goids_fin, "read_goids")
            return {'goids' : self.goids_fin}
        else:
            # Convert dict into 2D list retaining original section order
            sections_2d = []
            for section_name in self.sections_seen:
                if section_name in self.section2goids:
                    goids = self.section2goids.get(section_name)
                    chk_goids(goids, "GO IDs IN SECTION({S})".format(S=section_name))
                    sections_2d.append((section_name, goids))
            return {'sections' : sections_2d}

    def _prt_read_msg(self, prt, fin_txt, exclude_ungrouped):
        """Print which file was read and the number of GO IDs found."""
        if self.sections_seen or exclude_ungrouped:
            # dat = Grouper.get_summary_data(self.section2goids.items(), HdrgosSections.secdflt)
            dat = SummarySec2dHdrGos().summarize_sec2hdrgos(self.section2goids.items())
            sys.stdout.write(Grouper.fmtsum.format(
                GO_DESC='hdr', SECs=len(dat['S']), GOs=len(dat['G']),
                UNGRP="N/A", undesc="unused", ACTION="READ: ", FILE=fin_txt))
        elif self.goids_fin:
            prt.write("  {G} GO IDs READ: {FIN}\n".format(G=len(self.goids_fin), FIN=fin_txt))


# Copyright (C) 2016-2018, DV Klopfenstein, All rights reserved.
