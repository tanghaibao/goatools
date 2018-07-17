"""Product gene lists with ASCII art sections and GO IDs for each gene product."""

import sys
import collections as cx
from goatools.grouper.aart_geneproducts_one import AArtGeneProductSetsOne

__copyright__ = "Copyright (C) 2016-2018, DV Klopfenstein, H Tang, All rights reserved."
__author__ = "DV Klopfenstein"


class AArtGeneProductSetsAll(object):
    """Product gene lists with ASCII art sections and GO IDs for each gene product."""

    keys_exp = set([
        'fmtgo',        # Formatting for GO terms in grouped GO list
        'fmtgo2',       # Formatting for GO terms listed under each gene
        'fmtgene',      # Format for printing genes
        'fmtgene2',     # Format for printing genes in detailed report
        'itemid2name',  # Gene ID to Symbol
        ])

    #             A-Z                    a-z                     0-9
    #             : ; < = > ? @          !"#$&"()+,-./           [ \ ] ^ _
    #             { | } ~
    all_chrints = tuple(range(65, 91)) + tuple(range(97, 123)) + tuple(range(48, 58)) + \
                  tuple(range(58, 65)) + tuple(range(33, 48))  + tuple(range(91, 96)) + \
                  tuple(range(123, 127))

    def __init__(self, grprdflt, hdrobj, **kws):
        self.grprdflt = grprdflt  # GrouperDflts
        self.hdrobj = hdrobj      # HdrgosSections
        self.kws = {k:v for k, v in kws.items()}
        assert len(hdrobj.sections) <= len(self.all_chrints)
        self.sec2chr = cx.OrderedDict(
            [(s, chr(i)) for i, (s, _) in zip(self.all_chrints, hdrobj.sections)] + \
            [(hdrobj.secdflt, chr(self.all_chrints[len(hdrobj.sections)]))])
        self._init_kws()


    def run(self, name, goea_nts, log):
        """Run gene product ASCII art."""
        objaart = AArtGeneProductSetsOne(name, goea_nts, self)
        if self.hdrobj.sections:
            return objaart.prt_report_grp1(log)
        else:
            return objaart.prt_report_grp0(log)


    def prt_mrks(self, name_marks_list, prt=sys.stdout):
        """Print summary of all GOEAs.
           Example:
             Key for GO sections:
             A immune
             B viral/bacteria
             C neuro
             D cell death
             E lipid
             F adhesion
             G cell cycle
             H chromosome
             I development
             J extracellular matrix
             K ion
             L localization
             M membrane
             N metabolic
             O phosphorylation
             P signaling
             Q stimulus
             R prolif_differ
             S Misc.

             ABCDEFGHIJKLMNOPQRS
             XX.X..XXX..X.XX.XXX transient_increase
             XX.XXX.....X.X.XXXX consistent_increase
             XXXXXX..XXXXXXXXXXX late_increase
             ..X.....X.XX.X....X consistent_decrease
             ..X.XX..X.XX.XXX.XX late_decrease
        """
        if not name_marks_list:
            return
        # prt.write("\nKey for GO sections:\n")
        # self.prt_section_key(prt)
        prt.write("\n{HDR}\n".format(HDR=self.str_hdr()))
        for name, mark in name_marks_list:
            if mark is not None:
                prt.write("{MRKS} {NAME}\n".format(MRKS="".join(mark), NAME=name))
        prt.write("\n")

    def prt_section_key(self, prt=sys.stdout):
        """Print the section name and its alias."""
        for section_name, letter in self.sec2chr.items():
            prt.write("{ABC} {SEC_NAME}\n".format(ABC=letter, SEC_NAME=section_name))

    def str_hdr(self):
        """Return a string representing the section headers: """
        return "".join([c for _, c in self.sec2chr.items()])

    def get_chr2idx(self):
        """Return a dict with the ASCII art character as key and its index as value."""
        return {chr(ascii_int):idx for idx, ascii_int in enumerate(self.all_chrints)}

    def _init_kws(self):
        """Fill default values for keyword args, if necessary."""
        # Return user-specified GO formatting, if specfied:
        if 'fmtgo' not in self.kws:
            self.kws['fmtgo'] = self.grprdflt.gosubdag.prt_attr['fmt'] + "\n"
        if 'fmtgo2' not in self.kws:
            self.kws['fmtgo2'] = self.grprdflt.gosubdag.prt_attr['fmt'] + "\n"
        if 'fmtgene' not in self.kws:
            if 'itemid2name' not in self.kws:
                self.kws['fmtgene'] = "{AART} {ID}\n"
            else:
                self.kws['fmtgene'] = "{AART} {ID} {NAME}\n"
        if 'fmtgene2' not in self.kws:
            self.kws['fmtgene2'] = self.kws['fmtgene']

# Copyright (C) 2016-2018, DV Klopfenstein, H Tang, All rights reserved.
