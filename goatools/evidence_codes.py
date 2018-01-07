"""Manage evidence codes as reported by the Gene Ontology Consortium."""

__copyright__ = "Copyright (C) 2016-2018, DV Klopfenstein, H Tang. All rights reserved."
__author__ = "DV Klopfenstein"

import sys
import collections as cx

class EvidenceCodes(object):
    """From http://geneontology.org/page/guide-go-evidence-codes"""

    ntobj = cx.namedtuple("NtCode", "group name")

    code2name = cx.OrderedDict([
        # Experimental Evidence codes:
        ("EXP", ntobj._make(["Experimental", "Inferred from Experiment"])),
        ("IDA", ntobj._make(["Experimental", "Inferred from Direct Assay"])),
        ("IPI", ntobj._make(["Experimental", "Inferred from Physical Interaction"])),
        ("IMP", ntobj._make(["Experimental", "Inferred from Mutant Phenotype"])),
        ("IGI", ntobj._make(["Experimental", "Inferred from Genetic Interaction"])),
        ("IEP", ntobj._make(["Experimental", "Inferred from Expression Pattern"])),

        # Computational Analysis evidence codes
        ("ISS", ntobj._make(["Computational", "Inferred from Sequence or structural Similarity"])),
        ("ISO", ntobj._make(["Computational", "Inferred from Sequence Orthology"])),
        ("ISA", ntobj._make(["Computational", "Inferred from Sequence Alignment"])),
        ("ISM", ntobj._make(["Computational", "Inferred from Sequence Model"])),
        ("IGC", ntobj._make(["Computational", "Inferred from Genomic Context"])),
        ("IBA", ntobj._make(["Computational", "Inferred from Biological aspect of Ancestor"])),
        ("IBD", ntobj._make(["Computational", "Inferred from Biological aspect of Descendant"])),
        ("IKR", ntobj._make(["Computational", "Inferred from Key Residues"])),
        ("IRD", ntobj._make(["Computational", "Inferred from Rapid Divergence"])),
        ("RCA", ntobj._make(["Computational", "Inferred from Reviewed Computational Analysis"])),

        # Author Statement evidence codes
        ("TAS", ntobj._make(["Author", "Traceable Author Statement"])),
        ("NAS", ntobj._make(["Author", "Non-traceable Author Statement"])),

        # Curatorial Statement codes
        ("IC", ntobj._make(["Curatorial", "Inferred by Curator"])),
        ("ND", ntobj._make(["Curatorial", "No biological Data available"])),

        # Automatically-Assigned evidence code
        ("IEA", ntobj._make(["Automatic", "Inferred from Electronic Annotation"]))])

    def __init__(self):
        self.ev2idx = {ev:i for i, ev in enumerate(self.code2name.keys())}

    def sort_nts(self, nt_list, codekey):
        """Sort list of namedtuples such so evidence codes in same order as code2name."""
        # Problem is that some members in the nt_list do NOT have
        # codekey=EvidenceCode, then it returns None, which breaks py34 and 35
        # The fix here is that for these members, default to -1 (is this valid?)
        sortby = lambda nt: self.ev2idx.get(getattr(nt, codekey), -1)
        return sorted(nt_list, key=sortby)

    def get_grp_name(self, code):
        """Return group and name for an evidence code."""
        nt_code = self.code2name.get(code, None)
        if nt_code is not None:
            return nt_code.group, nt_code.name
        return "", ""

    def prt_ev_cnts(self, ctr, prt=sys.stdout):
        """Prints evidence code counts stored in a collections Counter."""
        for key, cnt in ctr.most_common():
            grp, name = self.get_grp_name(key.replace("NOT ", ""))
            prt.write("{CNT:7,} {EV:>7} {GROUP:<13} {NAME}\n".format(
                CNT=cnt, EV=key, GROUP=grp, NAME=name))

    def get_order(self, codes):
        """Return evidence codes in order shown in cod2name."""
        return sorted(codes, key=lambda e: [self.ev2idx.get(e)])

# Copyright (C) 2016-2018, DV Klopfenstein, H Tang. All rights reserved."
