"""Manages data from a DAVID Chart file."""

__copyright__ = "Copyright (C) 2016-2018, DV Klopfenstein, H Tang, All rights reserved."
__author__ = "DV Klopfenstein"

import os
import sys
import collections as cx
from goatools.wr_tbl import wr_xlsx
from goatools.wr_tbl import prt_txt

class DavidChartReader(object):
    """Manages data from a DAVID Chart file."""

    prt_flds = [
        'Category',
        'GO',
        'name',
        'Count',
        'Perc',
        'PValue',
        'List_Total',
        'Pop_Hits',
        'Pop_Total',
        'Fold_Enrichment',
        'Bonferroni',
        'Benjamini',
        'FDR',
        'Genes',
    ]

    fld2col_widths = {
        'Category':       10,  #  0         GOTERM_BP_ALL
        'GO':             11,  #  1 Term
        'name':           30,  #  1 Term
        'Count':           4,  #  2
        'Perc':            4,  #  3    %
        'PValue':          8,  #  4
        'Genes':          15,  #  5
        'List_Total':      4,  #  6
        'Pop_Hits':        4,  #  7
        'Pop_Total':       5,  #  8
        'Fold_Enrichment': 6,  #  9
        'Bonferroni':      8,  # 10
        'Benjamini':       8,  # 11
        'FDR':             8,  # 12
    }

    fld2prtfmt = {
        'Category':'{Category:10}',
        'GO':'{GO:10}',
        'name':'{name:30}',
        'Count':'{Count:8}',
        'Perc':'{Perc:4.1f}',
        'PValue':'{PValue:8.2e}',
        'Genes':'{Genes:30}',
        'List_Total':'{List_Total:8}',
        'Pop_Hits':'{Pop_Hits:8}',
        'Pop_Total':'{Pop_Total:8}',
        'Fold_Enrichment':'{Fold_Enrichment:6.4f}',
        'Bonferroni':'|{Bonferroni:8.2e}',
        'Benjamini':'{Benjamini:8.2e}',
        'FDR':'{FDR:8.2e}',
    }

    def __init__(self, fin_davidchart):
        self.fin_davidchart = fin_davidchart
        self.nts = _Init().get_nts(fin_davidchart)

    def wr_xlsx(self, fout_xlsx, nts):
        """Write specified namedtuples into an Excel spreadsheet."""
        wr_xlsx(fout_xlsx, nts, prt_flds=self.prt_flds, fld2col_widths=self.fld2col_widths)

    def prt_mdtbl(self, nts, prt_flds, prt=sys.stdout):
        """Write specified namedtuples into an Excel spreadsheet."""
        if prt_flds is None:
            prt_flds = self.prt_flds
        prtfmt = "|{DATA}|\n".format(DATA="|".join(self.fld2prtfmt[f] for f in prt_flds))
        prt_txt(prt, nts, prtfmt=prtfmt)

    def prt_num_sig(self, prt=sys.stdout, alpha=0.05):
        """Print the number of significant GO terms."""
        ctr = self.get_num_sig(alpha)
        prt.write("{N:6,} TOTAL: {TXT}\n".format(N=len(self.nts), TXT=" ".join([
            "FDR({FDR:4})".format(FDR=ctr['FDR']),
            "Bonferroni({B:4})".format(B=ctr['Bonferroni']),
            "Benjamini({B:4})".format(B=ctr['Benjamini']),
            "PValue({P:4})".format(P=ctr['PValue']),
            os.path.basename(self.fin_davidchart)])))

    def get_num_sig(self, alpha=0.05):
        """Print the number of significant results using various metrics."""
        # Get the number of significant GO terms
        ctr = cx.Counter()
        flds = set(['FDR', 'Bonferroni', 'Benjamini', 'PValue'])
        for ntd in self.nts:
            for fld in flds:
                if getattr(ntd, fld) < alpha:
                    ctr[fld] += 1
        return ctr


# pylint: disable=too-few-public-methods
class _Init(object):
    """Read DAVID Chart data and store in namedtuples."""

    nt_fields = [
        'Category',
        'GO',
        'name',
        'Count',
        'Perc',
        'PValue',
        'Genes',
        'Genes_set',
        'List_Total',
        'Pop_Hits',
        'Pop_Total',
        'Fold_Enrichment',
        'Bonferroni',
        'Benjamini',
        'FDR']

    def __init__(self):
        self.ntobj = cx.namedtuple("ntdavidchart", " ".join(self.nt_fields))

    def get_nts(self, fin_davidchart):
        """Read DAVID Chart file. Store each line in a namedtuple."""
        nts = []
        with open(fin_davidchart) as ifstrm:
            hdr_seen = False
            for line in ifstrm:
                line = line.rstrip()
                flds = line.split('\t')
                if hdr_seen:
                    ntd = self._init_nt(flds)
                    nts.append(ntd)
                else:
                    if line[:8] == 'Category':
                        assert len(flds) == 13, len(flds)
                        hdr_seen = True

            sys.stdout.write("  READ {N:5} GO IDs from DAVID Chart: {TSV}\n".format(
                N=len(nts), TSV=fin_davidchart))
        return nts

    def _init_nt(self, flds):
        """Given string fields from a DAVID chart file, return namedtuple."""
        vals = []
        # floats: %, PValue, Fold_Enrichment, Bonferroni, Benjamini, FDR
        floats = set([3, 4, 9, 10, 11, 12])
        # ints: Count, List_Total, Pop_Hits, Pop_Total
        ints = set([2, 6, 7, 8])
        for idx, valstr in enumerate(flds):
            if idx in floats:
                vals.append(float(valstr))
            elif idx in ints:
                vals.append(int(valstr))
            elif idx == 5:  # Genes
                vals.append(valstr)
                gene_set = valstr.split(', ')
                if gene_set and gene_set[0].isdigit():
                    gene_set = set(int(g) for g in gene_set)
                vals.append(gene_set)
            # Split 'Term' into GO and name
            elif idx == 1:
                goid, go_name = valstr.split('~')
                assert goid[:3] == "GO:" and len(goid) == 10
                vals.append(goid)
                vals.append(go_name)
            else:
                vals.append(valstr)
        return self.ntobj._make(vals)


# Copyright (C) 2016-2018, DV Klopfenstein, H Tang, All rights reserved.
