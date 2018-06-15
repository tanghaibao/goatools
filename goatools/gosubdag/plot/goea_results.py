"""Manages GO Term fill colors and bordercolors."""

__copyright__ = "Copyright (C) 2016-2017, DV Klopfenstein, H Tang, All rights reserved."
__author__ = "DV Klopfenstein"

import sys
import collections as cx


class GoeaResults(object):
    """Manages GOEA Results for plotting."""

    kws_set = set(['id2symbol', 'study_items', 'items_p_line ', 'pval_name'])

    dflt_items_p_line = 5 # study items (e.g., genes) per line on GO Terms
    fmtres = "{study_count} genes"

    alpha2col = cx.OrderedDict([
        # Enriched GOEA GO terms that are significant
        (0.005, 'mistyrose'),
        (0.010, 'moccasin'),
        (0.050, 'lemonchiffon1'),
        # GOEA GO terms that are not significant
        (1.000, 'grey95'),
    ])

    def __init__(self, goea_results, **kws):
        # kws: goea_results or go2nt
        assert goea_results, "NO GOEA RESULTS IN GoeaResults INPUTS"
        # GOATOOLs results as objects (WAS: Kws goea_results go2nt)
        self.go2res = {r.GO: r for r in goea_results}
        self.is_goterm = hasattr(goea_results[0], "_fldsdefprt")
        # GOATOOLs results as a list of namedtuples
        self.pval_name = self._init_pval_name(**kws)
        self.study_items = kws.get('study_items', None)
        self.study_items_max = self._init_study_items_max()
        self.items_p_line = kws['items_p_line'] if 'items_p_line' in kws else self.dflt_items_p_line
        self.id2symbol = kws['id2symbol'] if 'id2symbol' in kws else {}

    def prt_summary(self, prt=sys.stdout):
        """Print summary of GOEA plotting object."""
        desc = "NtGoeaResults" if self.is_goterm else "namedtuple"
        prt.write("{N} GOEA results from {O}. P-values stored in {P}.\n".format(
            N=len(self.go2res), O=desc, P=self.pval_name))

    def get_study_txt(self, goid):
        """Get GO text from GOEA study."""
        if goid in self.go2res:
            res = self.go2res[goid]
            if res.study_items is not None:
                return self._get_item_str(res)
            else:
                return self.fmtres.format(study_count=res.study_count)

    def set_goid2color_pval(self, goid2color):
        """Fill missing colors based on p-value of an enriched GO term."""
        alpha2col = self.alpha2col
        if self.pval_name is not None:
            pval_name = self.pval_name
            for goid, res in self.go2res.items():
                pval = getattr(res, pval_name, None)
                if pval is not None:
                    for alpha, color in alpha2col.items():
                        if pval <= alpha and res.study_count != 0:
                            if goid not in goid2color:
                                goid2color[goid] = color

    def get_goid2color_pval(self):
        """Return a go2color dict containing GO colors determined by P-value."""
        go2color = {}
        self.set_goid2color_pval(go2color)
        color_dflt = self.alpha2col[1.000]
        for goid in self.go2res:
            if goid not in go2color:
                go2color[goid] = color_dflt
        return go2color

    def _get_item_str(self, res):
        """Return genes in any of these formats:
              1. 19264, 17319, 12520, 12043, 74131, 22163, 12575
              2. Ptprc, Mif, Cd81, Bcl2, Sash3, Tnfrsf4, Cdkn1a
              3. 7: Ptprc, Mif, Cd81, Bcl2, Sash3...
        """
        ipl = self.items_p_line
        prt_items = sorted([self._get_genestr(itemid) for itemid in res.study_items])
        prt_multiline = [prt_items[i:i+ipl] for i in range(0, len(prt_items), ipl)]
        num_items = len(prt_items)
        if self.study_items_max is None:
            genestr = "\n".join([", ".join(str(e) for e in sublist) for sublist in prt_multiline])
            return "{N}) {GENES}".format(N=num_items, GENES=genestr)
        else:
            if num_items <= self.study_items_max:
                gene_lines = [", ".join(str(e) for e in sublist) for sublist in prt_multiline]
                genestr = "\n".join(gene_lines)
                return genestr
            else:
                short_list = prt_items[:self.study_items_max]
                short_mult = [short_list[i:i+ipl] for i in range(0, len(short_list), ipl)]
                short_lines = [", ".join(str(e) for e in sublist) for sublist in short_mult]
                short_str = "\n".join(short_lines)
                return "".join(["{N} genes; ".format(N=num_items), short_str, "..."])

    def _get_genestr(self, itemid):
        """Given a geneid, return the string geneid or a gene symbol."""
        if itemid in self.id2symbol:
            symbol = self.id2symbol[itemid]
            if symbol is not None:
                return symbol
        if isinstance(itemid, int):
            return str(itemid)
        return itemid


    def _init_pval_name(self, **kws):
        """Initialize pvalue attribute name."""
        if 'pval_name' in kws:
            return kws['pval_name']
        # If go2res contains GO Terms
        if self.is_goterm:
            return "p_{M}".format(M=next(iter(self.go2res.values())).get_method_name())
        # If go2res contains GO namedtuples
        for fld in next(iter(self.go2res.values()))._fields:
            if fld[:2] == 'p_' and fld != 'p_uncorrected':
                return fld

    def _init_study_items_max(self):
        """User can limit the number of genes printed in a GO term."""
        if self.study_items is None:
            return None
        if self.study_items is True:
            return None
        if isinstance(self.study_items, int):
            return self.study_items
        return None

# Copyright (C) 2016-2017, DV Klopfenstein, H Tang, All rights reserved.
