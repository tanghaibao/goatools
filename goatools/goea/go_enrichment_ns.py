# -*- coding: UTF-8 -*-
"""Runs Fisher's exact test, as well as multiple corrections for each of: BP MF CC"""

__copyright__ = "Copyright (C) 2010-2019, H Tang et al., All rights reserved."
__author__ = "various"

import itertools
from goatools.go_enrichment import GOEnrichmentStudy


class GOEnrichmentStudyNS:
    """Runs Fisher's exact test, as well as multiple corrections for each of: BP MF CC"""

    # pylint: disable=too-many-arguments
    def __init__(self, pop, ns2assoc, godag, propagate_counts=True, alpha=.05, methods=None, **kws):
        self.ns2objgoea = self._ns2o(pop, ns2assoc, godag, propagate_counts, alpha, methods, **kws)

    def run_study(self, study_ids):
        """Run GOEAs for each namespace, BP MF CC"""
        ns2results = {ns:o.run_study(study_ids) for ns, o in sorted(self.ns2objgoea.items())}
        return list(itertools.chain.from_iterable(ns2results.values()))

    def wr_xlsx(self, fout_xlsx, goea_results, **kws):
        """Write to spreadsheet format"""
        next(iter(self.ns2objgoea.values())).wr_xlsx(fout_xlsx, goea_results, **kws)

    def wr_tsv(self, fout_tsv, goea_results, **kws):
        """Write to spreadsheet format"""
        next(iter(self.ns2objgoea.values())).wr_tsv(fout_tsv, goea_results, **kws)

    def wr_txt(self, fout_tsv, goea_results, **kws):
        """Write to spreadsheet format"""
        next(iter(self.ns2objgoea.values())).wr_txt(fout_tsv, goea_results, **kws)

    @staticmethod
    def _ns2o(pop, ns2assoc, godag, propagate_counts, alpha, methods, **kws):
        return {
            ns:GOEnrichmentStudy(pop, a, godag, propagate_counts, alpha, methods, name=ns, **kws) \
                for ns, a in sorted(ns2assoc.items())}

    def get_assoc(self):
        """Get a list of all GO sets in all (BP, MF, CC) associations"""
        return {geneid:gos for o in self.ns2objgoea.values() for geneid, gos in o.assoc.items()}


# Copyright (C) 2010-2019, H Tang et al., All rights reserved.
