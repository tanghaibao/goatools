"""Manage evidence codes as reported by the Gene Ontology Consortium."""

__copyright__ = "Copyright (C) 2016, DV Klopfenstein, H Tang. All rights reserved."
__author__ = "DV Klopfenstein"

import sys
import collections as cx

class EvidenceCodes(object):
    """From http://geneontology.org/page/guide-go-evidence-codes"""

    ntobj = cx.namedtuple("NtCode", "group name")

    code2name = {
        # Experimental Evidence codes:
        "EXP":ntobj(group="Experimental", name="Inferred from Experiment"),
        "IDA":ntobj(group="Experimental", name="Inferred from Direct Assay"),
        "IPI":ntobj(group="Experimental", name="Inferred from Physical Interaction"),
        "IMP":ntobj(group="Experimental", name="Inferred from Mutant Phenotype"),
        "IGI":ntobj(group="Experimental", name="Inferred from Genetic Interaction"),
        "IEP":ntobj(group="Experimental", name="Inferred from Expression Pattern"),

        # Computational Analysis evidence codes:
        "ISS":ntobj(group="Computational", name="Inferred from Sequence or structural Similarity"),
        "ISO":ntobj(group="Computational", name="Inferred from Sequence Orthology"),
        "ISA":ntobj(group="Computational", name="Inferred from Sequence Alignment"),
        "ISM":ntobj(group="Computational", name="Inferred from Sequence Model"),
        "IGC":ntobj(group="Computational", name="Inferred from Genomic Context"),
        "IBA":ntobj(group="Computational", name="Inferred from Biological aspect of Ancestor"),
        "IBD":ntobj(group="Computational", name="Inferred from Biological aspect of Descendant"),
        "IKR":ntobj(group="Computational", name="Inferred from Key Residues"),
        "IRD":ntobj(group="Computational", name="Inferred from Rapid Divergence"),
        "RCA":ntobj(group="Computational", name="Inferred from Reviewed Computational Analysis"),

        # Author Statement evidence codes:
        "TAS":ntobj(group="Author", name="Traceable Author Statement"),
        "NAS":ntobj(group="Author", name="Non-traceable Author Statement"),

        # Curatorial Statement codes:
        "IC":ntobj(group="Curatorial", name="Inferred by Curator"),
        "ND":ntobj(group="Curatorial", name="No biological Data available"),

        # Automatically-Assigned evidence code:
        "IEA":ntobj(group="Automatic", name="Inferred from Electronic Annotation")}

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

# Copyright (C) 2016, DV Klopfenstein, H Tang. All rights reserved."
