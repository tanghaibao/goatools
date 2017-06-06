#!/usr/bin/env python
# -*- coding: UTF-8 -*-

"""
python %prog study.file population.file gene-association.file

This program returns P-values for functional enrichment in a cluster of
study genes using Fisher's exact test, and corrected for multiple testing
(including Bonferroni, Holm, Sidak, and false discovery rate)
"""

from __future__ import absolute_import

__copyright__ = "Copyright (C) 2010-2017, H Tang et al., All rights reserved."
__author__ = "various"

import sys
import collections as cx
import datetime
from functools import partial
import multiprocessing


from goatools.multiple_testing import Methods, Bonferroni, Sidak, HolmBonferroni, FDR, calc_qval
from goatools.ratio import get_terms, count_terms, is_ratio_different
import goatools.wr_tbl as RPT
from goatools.pvalcalc import FisherFactory




def compute_pvals(allterms, calc_pvalue, go2studyitems, go2popitems,
                  study_n, pop_n):

    import os
    print("Pid=", os.getpid())

    results = {}
    for term in allterms:

        study_items = go2studyitems.get(term, set())
        study_count = len(study_items)
        pop_items = go2popitems.get(term, set())
        pop_count = len(pop_items)

        results[term] = calc_pvalue(study_count, study_n, pop_count, pop_n)

    return results


class GOEnrichmentRecord(object):
    """Represents one result (from a single GOTerm) in the GOEnrichmentStudy
    """
    namespace2NS = cx.OrderedDict([
        ('biological_process', 'BP'),
        ('molecular_function', 'MF'),
        ('cellular_component', 'CC')])

    # Fields seen in every enrichment result
    _fldsdefprt = [
        "GO",
        "NS",
        "enrichment",
        "name",
        "ratio_in_study",
        "ratio_in_pop",
        "p_uncorrected",
        "depth",
        "study_count",
        "study_items"]
    _fldsdeffmt = ["%s"]*3 + ["%-30s"] + ["%d/%d"] * 2 + ["%.3g"] + ["%d"] * 2 + ["%15s"]

    _flds = set(_fldsdefprt).intersection(
        set(['study_items', 'study_count', 'study_n', 'pop_items', 'pop_count', 'pop_n']))

    def __init__(self, **kwargs):
        # Methods seen in current enrichment result
        self._methods = []
        for k, v in kwargs.items():
            setattr(self, k, v)
            if k == 'ratio_in_study':
                setattr(self, 'study_count', v[0])
                setattr(self, 'study_n', v[1])
            if k == 'ratio_in_pop':
                setattr(self, 'pop_count', v[0])
                setattr(self, 'pop_n', v[1])
        self._init_enrichment()
        self.goterm = None  # the reference to the GOTerm

    def get_method_name(self):
        """Return name of first method in the _methods list."""
        return self._methods[0].fieldname

    def get_pvalue(self):
        """Returns pval for 1st method, if it exists. Else returns uncorrected pval."""
        if self._methods:
            return getattr(self, "p_{m}".format(m=self.get_method_name()))
        return getattr(self, "p_uncorrected")

    def set_corrected_pval(self, nt_method, pvalue):
        """Add object attribute based on method name."""
        self._methods.append(nt_method)
        fieldname = "".join(["p_", nt_method.fieldname])
        setattr(self, fieldname, pvalue)

    def __str__(self, indent=False):
        field_data = [getattr(self, f, "n.a.") for f in self._fldsdefprt[:-1]] + \
                     [getattr(self, "p_{}".format(m.fieldname)) for m in self._methods] + \
                     [", ".join(sorted(getattr(self, self._fldsdefprt[-1], set())))]
        fldsdeffmt = self._fldsdeffmt
        field_formatter = fldsdeffmt[:-1] + ["%.3g"]*len(self._methods) + [fldsdeffmt[-1]]
        self._chk_fields(field_data, field_formatter)

        # default formatting only works for non-"n.a" data
        for i, f in enumerate(field_data):
            if f == "n.a.":
                field_formatter[i] = "%s"

        # print dots to show the level of the term
        dots = self.get_indent_dots() if indent else ""
        prtdata = "\t".join(a % b for (a, b) in zip(field_formatter, field_data))
        return "".join([dots, prtdata])

    def get_indent_dots(self):
        """Get a string of dots ("....") representing the level of the GO term."""
        return "." * self.goterm.level if self.goterm is not None else ""

    @staticmethod
    def _chk_fields(field_data, field_formatter):
        """Check that expected fields are present."""
        if len(field_data) == len(field_formatter):
            return
        len_dat = len(field_data)
        len_fmt = len(field_formatter)
        msg = [
            "FIELD DATA({d}) != FORMATTER({f})".format(d=len_dat, f=len_fmt),
            "DAT({N}): {D}".format(N=len_dat, D=field_data),
            "FMT({N}): {F}".format(N=len_fmt, F=field_formatter)]
        raise Exception("\n".join(msg))

    def __repr__(self):
        return "GOEnrichmentRecord({GO})".format(GO=self.GO)

    def set_goterm(self, goid):
        """Set goterm and copy GOTerm's name and namespace."""
        self.goterm = goid.get(self.GO, None)
        present = self.goterm is not None
        self.name = self.goterm.name if present else "n.a."
        self.NS = self.namespace2NS[self.goterm.namespace] if present else "XX"

    def _init_enrichment(self):
        """Mark as 'enriched' or 'purified'."""
        self.enrichment = 'e' if ((1.0 * self.study_count / self.study_n) >
                                  (1.0 * self.pop_count / self.pop_n)) else 'p'

    def update_remaining_fldsdefprt(self, min_ratio=None):
        """Finish updating self (GOEnrichmentRecord) field, is_ratio_different."""
        self.is_ratio_different = is_ratio_different(min_ratio, self.study_count,
                                                     self.study_n, self.pop_count, self.pop_n)


    # -------------------------------------------------------------------------------------
    # Methods for getting flat namedtuple values from GOEnrichmentRecord object
    def get_prtflds_default(self):
        """Get default fields."""
        return self._fldsdefprt[:-1] + \
               ["p_{M}".format(M=m.fieldname) for m in self._methods] + \
               [self._fldsdefprt[-1]]

    def get_prtflds_all(self):
        """When converting to a namedtuple, get all possible fields in their original order."""
        flds = []
        dont_add = set(['_parents', '_methods'])
        # Fields: GO NS enrichment name ratio_in_study ratio_in_pop p_uncorrected
        #         depth study_count p_sm_bonferroni p_fdr_bh study_items
        self._flds_append(flds, self.get_prtflds_default(), dont_add)
        # Fields: GO NS goterm
        #         ratio_in_pop pop_n pop_count pop_items name
        #         ratio_in_study study_n study_count study_items
        #         _methods enrichment p_uncorrected p_sm_bonferroni p_fdr_bh
        self._flds_append(flds, vars(self).keys(), dont_add)
        # Fields: name level is_obsolete namespace id depth parents children _parents alt_ids
        self._flds_append(flds, vars(self.goterm).keys(), dont_add)
        return flds

    @staticmethod
    def _flds_append(flds, addthese, dont_add):
        """Retain order of fields as we add them once to the list."""
        for fld in addthese:
            if fld not in flds and fld not in dont_add:
                flds.append(fld)

    def get_field_values(self, fldnames, rpt_fmt=True):
        """Get flat namedtuple fields for one GOEnrichmentRecord."""
        row = []
        # Loop through each user field desired
        for fld in fldnames:
            # 1. Check the GOEnrichmentRecord's attributes
            val = getattr(self, fld, None)
            if val is not None:
                if rpt_fmt:
                    val = self._get_rpt_fmt(fld, val)
                row.append(val)
            else:
                # 2. Check the GO object for the field
                val = getattr(self.goterm, fld, None)
                if rpt_fmt:
                    val = self._get_rpt_fmt(fld, val)
                if val is not None:
                    row.append(val)
                else:
                    # 3. Field not found, raise Exception
                    self._err_fld(fld, fldnames, row)
            if rpt_fmt:
                assert not isinstance(val, list), \
                   "UNEXPECTED LIST: FIELD({F}) VALUE({V}) FMT({P})".format(
                       P=rpt_fmt, F=fld, V=val)
        return row

    @staticmethod
    def _get_rpt_fmt(fld, val):
        """Return values in a format amenable to printing in a table."""
        if fld.startswith("ratio_"):
            return "{N}/{TOT}".format(N=val[0], TOT=val[1])
        elif fld in set(['study_items', 'pop_items', 'alt_ids']):
            return ", ".join([str(v) for v in sorted(val)])
        return val

    def _err_fld(self, fld, fldnames):
        """Unrecognized field. Print detailed Failure message."""
        msg = ['ERROR. UNRECOGNIZED FIELD({F})'.format(F=fld)]
        actual_flds = set(self.get_prtflds_default() + self.goterm.__dict__.keys())
        bad_flds = set(fldnames).difference(set(actual_flds))
        if bad_flds:
            msg.append("\nGOEA RESULT FIELDS: {}".format(" ".join(self._fldsdefprt)))
            msg.append("GO FIELDS: {}".format(" ".join(self.goterm.__dict__.keys())))
            msg.append("\nFATAL: {N} UNEXPECTED FIELDS({F})\n".format(
                N=len(bad_flds), F=" ".join(bad_flds)))
            msg.append("  {N} User-provided fields:".format(N=len(fldnames)))
            for idx, fld in enumerate(fldnames, 1):
                mrk = "ERROR -->" if fld in bad_flds else ""
                msg.append("  {M:>9} {I:>2}) {F}".format(M=mrk, I=idx, F=fld))
        raise Exception("\n".join(msg))


class GOEnrichmentStudy(object):
    """Runs Fisher's exact test, as well as multiple corrections
    """
    # Default Excel table column widths for GOEA results
    default_fld2col_widths = {
        'NS'        :  3,
        'GO'        : 12,
        'level'     :  3,
        'enrichment':  1,
        'name'      : 60,
        'ratio_in_study':  8,
        'ratio_in_pop'  : 12,
        'study_items'   : 15,
    }

    def __init__(self, pop, assoc, obo_dag, propagate_counts=True, alpha=.05, methods=None, **kws):
        self.log = kws['log'] if 'log' in kws else sys.stdout
        self._run_multitest = {
            'local':lambda iargs: self._run_multitest_local(iargs),
            'statsmodels':lambda iargs: self._run_multitest_statsmodels(iargs)}
        self.pop = pop
        self.pop_n = len(pop)
        self.assoc = assoc
        self.obo_dag = obo_dag
        self.alpha = alpha
        if methods is None:
            methods = ["bonferroni", "sidak", "holm"]
        self.methods = Methods(methods)
        self.pval_obj = FisherFactory(**kws).pval_obj

        if propagate_counts:
            sys.stderr.write("Propagating term counts to parents ..\n")
            obo_dag.update_association(assoc)
        self.go2popitems = get_terms("population", pop, assoc, obo_dag, self.log)

    def run_study(self, study, **kws):
        """Run Gene Ontology Enrichment Study (GOEA) on study ids."""
        # Key-word arguments:
        methods = Methods(kws['methods']) if 'methods' in kws else self.methods
        alpha = kws['alpha'] if 'alpha' in kws else self.alpha
        log = kws['log'] if 'log' in kws else self.log
        # Calculate uncorrected pvalues
        results = self._get_pval_uncorr(study, log)
        if not results:
            return []

        # Do multipletest corrections on uncorrected pvalues and update results
        self._run_multitest_corr(results, methods, alpha, study, log)

        for rec in results:
            # get go term for name and level
            rec.set_goterm(self.obo_dag)

        # 'keep_if' can be used to keep only significant GO terms. Example:
        #     >>> keep_if = lambda nt: nt.p_fdr_bh < 0.05 # if results are significant
        #     >>> goea_results = goeaobj.run_study(geneids_study, keep_if=keep_if)
        if 'keep_if' in kws:
            keep_if = kws['keep_if']
            results = [r for r in results if keep_if(r)]

        # Default sort order: First, sort by BP, MF, CC. Second, sort by pval
        results.sort(key=lambda r: [r.NS, r.p_uncorrected])

        if log is not None:
            log.write("  {MSG}\n".format(MSG="\n  ".join(self.get_results_msg(results, study))))

        return results # list of GOEnrichmentRecord objects

    def run_study_nts(self, study, **kws):
        """Run GOEA on study ids. Return results as a list of namedtuples."""
        goea_results = self.run_study(study, **kws)
        return get_goea_nts_all(goea_results)

    def get_results_msg(self, results, study):
        """Return summary for GOEA results."""
        # To convert msg list to string: "\n".join(msg)
        msg = []
        if results:
            stu_items, num_gos_stu = self.get_item_cnt(results, "study_items")
            pop_items, num_gos_pop = self.get_item_cnt(results, "pop_items")
            msg.append("{M:,} GO terms are associated with {N:,} of {NT:,} study items".format(
                N=len(stu_items), NT=len(set(study)), M=num_gos_stu))
            msg.append("{M:,} GO terms are associated with {N:,} of {NT:,} population items".format(
                N=len(pop_items), NT=self.pop_n, M=num_gos_pop))
        return msg

    def _get_pval_uncorr(self, study, log=sys.stdout):
        """Calculate the uncorrected pvalues for study items."""
        if log is not None:
            log.write("Calculating uncorrected p-values using {PFNC}\n".format(PFNC=self.pval_obj.name))
        go2studyitems = get_terms("study", study, self.assoc, self.obo_dag, log)
        pop_n, study_n = self.pop_n, len(study)
        allterms = set(go2studyitems.keys()).union(
            set(self.go2popitems.keys()))

        # if self.pval_obj.log is a file handle, which we can not serialize, so we could
        # not transfer self.pval_obj.calc_pvalue to another python process with multiprocessing.
        # there fore we "path" the object which will later be restored again.
        old = self.pval_obj.log
        self.pval_obj.log = None
        calc_pvalue = self.pval_obj.calc_pvalue

        # -1 avoids freezing of the machine:
        n_procs = multiprocessing.cpu_count() - 1

        p = multiprocessing.Pool(n_procs)
        n = len(allterms)

        allterms = list(allterms)
        fragments = [allterms[i::n_procs] for i in range(n_procs)]

        remote_func = partial(compute_pvals, calc_pvalue=calc_pvalue, go2studyitems=go2studyitems,
                              go2popitems=self.go2popitems, study_n=study_n, pop_n=pop_n)

        all_p_values = p.map(remote_func, fragments)

        # restore patched file handle
        self.pval_obj.log = old

        results = []

        for p_values in all_p_values:

            for term, p_value in p_values.items():

                study_items = go2studyitems.get(term, set())
                study_count = len(study_items)
                pop_items = self.go2popitems.get(term, set())
                pop_count = len(pop_items)

                one_record = GOEnrichmentRecord(
                    GO=term,
                    p_uncorrected=p_value,
                    study_items=study_items,
                    pop_items=pop_items,
                    ratio_in_study=(study_count, study_n),
                    ratio_in_pop=(pop_count, pop_n))

                results.append(one_record)

        return results

    def _run_multitest_corr(self, results, usr_methods, alpha, study, log):
        """Do multiple-test corrections on uncorrected pvalues."""
        assert 0 < alpha < 1, "Test-wise alpha must fall between (0, 1)"
        pvals = [r.p_uncorrected for r in results]
        NtMt = cx.namedtuple("NtMt", "results pvals alpha nt_method study")

        for nt_method in usr_methods:
            ntmt = NtMt(results, pvals, alpha, nt_method, study)
            if log is not None:
                log.write("Running multitest correction: {MSRC} {METHOD}\n".format(
                    MSRC=ntmt.nt_method.source, METHOD=ntmt.nt_method.method))
            self._run_multitest[nt_method.source](ntmt)

    def _run_multitest_statsmodels(self, ntmt):
        """Use multitest mthods that have been implemented in statsmodels."""
        # Only load statsmodels if it is used
        multipletests = self.methods.get_statsmodels_multipletests()
        results = multipletests(ntmt.pvals, ntmt.alpha, ntmt.nt_method.method)
        pvals_corrected = results[1] # reject_lst, pvals_corrected, alphacSidak, alphacBonf
        self._update_pvalcorr(ntmt, pvals_corrected)

    def _run_multitest_local(self, ntmt):
        """Use multitest mthods that have been implemented locally."""
        corrected_pvals = None
        method = ntmt.nt_method.method
        if method == "bonferroni":
            corrected_pvals = Bonferroni(ntmt.pvals, ntmt.alpha).corrected_pvals
        elif method == "sidak":
            corrected_pvals = Sidak(ntmt.pvals, ntmt.alpha).corrected_pvals
        elif method == "holm":
            corrected_pvals = HolmBonferroni(ntmt.pvals, ntmt.alpha).corrected_pvals
        elif method == "fdr":
            # get the empirical p-value distributions for FDR
            term_pop = getattr(self, 'term_pop', None)
            if term_pop is None:
                term_pop = count_terms(self.pop, self.assoc, self.obo_dag)
            p_val_distribution = calc_qval(len(ntmt.study),
                                           self.pop_n,
                                           self.pop, self.assoc,
                                           term_pop, self.obo_dag)
            corrected_pvals = FDR(p_val_distribution,
                                  ntmt.results, ntmt.alpha).corrected_pvals

        self._update_pvalcorr(ntmt, corrected_pvals)

    @staticmethod
    def _update_pvalcorr(ntmt, corrected_pvals):
        """Add data members to store multiple test corrections."""
        if corrected_pvals is None:
            return
        for rec, val in zip(ntmt.results, corrected_pvals):
            rec.set_corrected_pval(ntmt.nt_method, val)

    # Methods for writing results into tables: text, tab-separated, Excel spreadsheets
    def wr_txt(self, fout_txt, goea_results, prtfmt=None, **kws):
        """Print GOEA results to text file."""
        if not goea_results:
            sys.stdout.write("      0 GOEA results. NOT WRITING {FOUT}\n".format(FOUT=fout_txt))
            return
        with open(fout_txt, 'w') as prt:
            data_nts = self.prt_txt(prt, goea_results, prtfmt, **kws)
            log = self.log if self.log is not None else sys.stdout
            log.write("  {N:>5} GOEA results for {CUR:5} study items. WROTE: {F}\n".format(
                N=len(data_nts),
                CUR=len(get_study_items(goea_results)),
                F=fout_txt))

    def prt_txt(self, prt, goea_results, prtfmt=None, **kws):
        """Print GOEA results in text format."""
        if prtfmt is None:
            prtfmt = "{GO} {NS} {p_uncorrected:5.2e} {study_count:>5} {name}\n"
        prtfmt = self.adjust_prtfmt(prtfmt)
        prt_flds = RPT.get_fmtflds(prtfmt)
        data_nts = get_goea_nts_prt(goea_results, prt_flds, **kws)
        RPT.prt_txt(prt, data_nts, prtfmt, prt_flds, **kws)
        return data_nts

    def wr_xlsx(self, fout_xlsx, goea_results, **kws):
        """Write a xlsx file."""
        # kws: prt_if indent
        prt_flds = kws.get('prt_flds', self.get_prtflds_default(goea_results))
        xlsx_data = get_goea_nts_prt(goea_results, prt_flds, **kws)
        if 'fld2col_widths' not in kws:
            kws['fld2col_widths'] = {f:self.default_fld2col_widths.get(f, 8) for f in prt_flds}
        RPT.wr_xlsx(fout_xlsx, xlsx_data, **kws)

    def wr_tsv(self, fout_tsv, goea_results, **kws):
        """Write tab-separated table data to file"""
        prt_flds = kws.get('prt_flds', self.get_prtflds_default(goea_results))
        tsv_data = get_goea_nts_prt(goea_results, prt_flds, **kws)
        RPT.wr_tsv(fout_tsv, tsv_data, **kws)

    def prt_tsv(self, prt, goea_results, **kws):
        """Write tab-separated table data"""
        prt_flds = kws.get('prt_flds', self.get_prtflds_default(goea_results))
        tsv_data = get_goea_nts_prt(goea_results, prt_flds, **kws)
        RPT.prt_tsv(prt, tsv_data, prt_flds, **kws)

    @staticmethod
    def adjust_prtfmt(prtfmt):
        """Adjust format_strings for legal values."""
        prtfmt = prtfmt.replace("{p_holm-sidak", "{p_holm_sidak")
        prtfmt = prtfmt.replace("{p_simes-hochberg", "{p_simes_hochberg")
        return prtfmt

    @staticmethod
    def get_NS2nts(results, fldnames=None, **kws):
        """Get namedtuples of GOEA results, split into BP, MF, CC."""
        NS2nts = cx.defaultdict(list)
        nts = get_goea_nts_all(results, fldnames, **kws)
        for nt in nts:
            NS2nts[nt.NS].append(nt)
        return NS2nts

    @staticmethod
    def get_item_cnt(results, attrname="study_items"):
        """Get all study or population items (e.g., geneids)."""
        items = set()
        go_cnt = 0
        for rec in results:
            if hasattr(rec, attrname):
                items_cur = getattr(rec, attrname)
                # Only count GO term if there are items in the set.
                if len(items_cur) != 0:
                    items |= items_cur
                    go_cnt += 1
        return items, go_cnt

    @staticmethod
    def get_prtflds_default(results):
        """Get default fields names. Used in printing GOEA results.

           Researchers can control which fields they want to print in the GOEA results
           or they can use the default fields.
        """
        if results:
            return results[0].get_prtflds_default()
        return []

    @staticmethod
    def print_summary(results, min_ratio=None, indent=False, pval=0.05):
        """Print summary."""
        from .version import __version__ as version

        # Header contains provenance and parameters
        print("# Generated by GOATOOLS v{0} ({1})".format(version, datetime.date.today()))
        print("# min_ratio={0} pval={1}".format(min_ratio, pval))

        # field names for output
        if results:
            print("\t".join(GOEnrichmentStudy.get_prtflds_default(results)))

        for rec in results:
            # calculate some additional statistics
            # (over_under, is_ratio_different)
            rec.update_remaining_fldsdefprt(min_ratio=min_ratio)

            if pval is not None and rec.p_uncorrected >= pval:
                continue

            if rec.is_ratio_different:
                print(rec.__str__(indent=indent))

    def wr_py_goea_results(self, fout_py, goea_results, **kws):
        """Save GOEA results into Python package containing list of namedtuples."""
        var_name = kws.get("var_name", "goea_results")
        docstring = kws.get("docstring", "")
        sortby = kws.get("sortby", None)
        if goea_results:
            from goatools.nt_utils import wr_py_nts
            nts_goea = goea_results
            # If list has GOEnrichmentRecords or verbose namedtuples, exclude some fields.
            if hasattr(goea_results[0], "_fldsdefprt") or hasattr(goea_results[0], 'goterm'):
                # Exclude some attributes from the namedtuple when saving results
                # to a Python file because the information is redundant or verbose.
                nts_goea = get_goea_nts_prt(goea_results)
            docstring = "\n".join([docstring, "# {VER}\n\n".format(VER=self.obo_dag.version)])
            assert hasattr(nts_goea[0], '_fields')
            if sortby is None:
                sortby = lambda nt: getattr(nt, 'p_uncorrected')
            nts_goea = sorted(nts_goea, key=sortby)
            wr_py_nts(fout_py, nts_goea, docstring, var_name)

def get_study_items(goea_results):
    """Get all study items (e.g., geneids)."""
    study_items = set()
    for rec in goea_results:
        study_items |= rec.study_items
    return study_items

def get_goea_nts_prt(goea_results, fldnames=None, **usr_kws):
    """Return list of namedtuples removing fields which are redundant or verbose."""
    kws = usr_kws.copy()
    if 'not_fldnames' not in kws:
        kws['not_fldnames'] = ['goterm', 'parents', 'children', 'id']
    if 'rpt_fmt' not in kws:
        kws['rpt_fmt'] = True
    return get_goea_nts_all(goea_results, fldnames, **kws)

def get_goea_nts_all(goea_results, fldnames=None, **kws):
    """Get namedtuples containing user-specified (or default) data from GOEA results.

        Reformats data from GOEnrichmentRecord objects into lists of
        namedtuples so the generic table writers may be used.
    """
    data_nts = [] # A list of namedtuples containing GOEA results
    if not goea_results:
        return data_nts
    keep_if = kws.get('keep_if', None)
    rpt_fmt = kws.get('rpt_fmt', False)
    indent = kws.get('indent', False)
    # I. FIELD (column) NAMES
    not_fldnames = kws.get('not_fldnames', None)
    if fldnames is None:
        fldnames = get_fieldnames(goea_results[0])
    # Ia. Explicitly exclude specific fields from named tuple
    if not_fldnames is not None:
        fldnames = [f for f in fldnames if f not in not_fldnames]
    nttyp = cx.namedtuple("NtGoeaResults", " ".join(fldnames))
    goid_idx = fldnames.index("GO") if 'GO' in fldnames else None
    # II. Loop through GOEA results stored in a GOEnrichmentRecord object
    for goerec in goea_results:
        vals = get_field_values(goerec, fldnames, rpt_fmt)
        if indent:
            vals[goid_idx] = "".join([goerec.get_indent_dots(), vals[goid_idx]])
        ntobj = nttyp._make(vals)
        if keep_if is None or keep_if(ntobj):
            data_nts.append(ntobj)
    return data_nts

def get_field_values(item, fldnames, rpt_fmt=None):
    """Return fieldnames and values of either a namedtuple or GOEnrichmentRecord."""
    if hasattr(item, "_fldsdefprt"): # Is a GOEnrichmentRecord
        return item.get_field_values(fldnames, rpt_fmt)
    if hasattr(item, "_fields"): # Is a namedtuple
        return [getattr(item, f) for f in fldnames]

def get_fieldnames(item):
    """Return fieldnames of either a namedtuple or GOEnrichmentRecord."""
    if hasattr(item, "_fldsdefprt"): # Is a GOEnrichmentRecord
        return item.get_prtflds_all()
    if hasattr(item, "_fields"): # Is a namedtuple
        return item._fields

# Copyright (C) 2010-2017, H Tang et al., All rights reserved.
