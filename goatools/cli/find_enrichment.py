# -*- coding: UTF-8 -*-
"""
python find_enrichment.py study.file population.file gene-association.file

This program returns P-values for functional enrichment in a cluster of study
genes using Fisher's exact test, and corrected for multiple testing (including
Bonferroni, Holm, Sidak, and false discovery rate).

About significance cutoff:
--alpha: test-wise alpha; for each GO term, what significance level to apply
        (most often you don't need to change this other than 0.05 or 0.01)
--pval: experiment-wise alpha; for the entire experiment, what significance
        level to apply after Bonferroni correction
"""

from __future__ import print_function

__copyright__ = "Copyright (C) 2010-2019, H Tang et al. All rights reserved."
__author__ = "various"

import os
import sys
import re
import argparse

from goatools.evidence_codes import EvidenceCodes

from goatools.obo_parser import GODag
from goatools.goea.go_enrichment_ns import GOEnrichmentStudyNS
from goatools.multiple_testing import Methods
from goatools.pvalcalc import FisherFactory
from goatools.rpt.goea_nt_xfrm import MgrNtGOEAs
from goatools.rpt.prtfmt import PrtFmt
from goatools.semantic import TermCounts
from goatools.wr_tbl import prt_tsv_sections
from goatools.godag.consts import RELATIONSHIP_SET
from goatools.godag.consts import chk_relationships
from goatools.godag.prtfncs import GoeaPrintFunctions
from goatools.anno.factory import get_anno_desc
from goatools.anno.factory import get_objanno
from goatools.cli.gos_get import GetGOs

from goatools.gosubdag.gosubdag import GoSubDag
from goatools.grouper.read_goids import read_sections
from goatools.grouper.grprdflts import GrouperDflts
from goatools.grouper.hdrgos import HdrgosSections
from goatools.grouper.grprobj import Grouper
from goatools.grouper.sorter import Sorter
from goatools.grouper.aart_geneproducts_all import AArtGeneProductSetsAll
from goatools.grouper.wr_sections import WrSectionsTxt
from goatools.grouper.wrxlsx import WrXlsxSortedGos
OBJPRTRES = GoeaPrintFunctions()


# pylint: disable=too-few-public-methods
class GoeaCliArgs:
    """Extracts arguments from the command-line."""

    def __init__(self):
        self.args = self._init_args()

    def _init_args(self):
        """Get enrichment arg parser."""

        #pylint: disable=invalid-name
        p = argparse.ArgumentParser(__doc__,
                                    formatter_class=argparse.ArgumentDefaultsHelpFormatter)

        p.add_argument('filenames', type=str, nargs=3,
                       help='data/study data/population data/association')
        p.add_argument('--annofmt', default=None, type=str,
                       help=('Annotation file format. '
                             'Not needed if type can be determined using filename'),
                       choices=['gene2go', 'gaf', 'gpad', 'id2gos'])
        p.add_argument('--taxid', default=9606, type=int,
                       help="When using NCBI's gene2go annotation file, specify desired taxid")
        p.add_argument('--alpha', default=0.05, type=float,
                       help='Test-wise alpha for multiple testing')
        p.add_argument('--pval', default=.05, type=float,
                       help='Only print results with uncorrected p-value < PVAL.')
        p.add_argument('--pval_field', type=str,
                       help='Only print results when PVAL_FIELD < PVAL.')
        p.add_argument('--outfile', default=None, type=str,
                       help='Write enrichment results into xlsx or tsv file')
        p.add_argument('--ns', default='BP,MF,CC', type=str,
                       help='Limit GOEA to specified branch categories. '
                            'BP=Biological Process; '
                            'MF=Molecular Function; '
                            'CC=Cellular Component')
        p.add_argument('--id2sym', default=None, type=str,
                       help='ASCII file containing one geneid and its symbol per line')
        p.add_argument('--sections', default=None, type=str,
                       help=('Use sections file for printing grouped GOEA results. '
                             'Example SECTIONS values:\n'
                             'goatools.test_data.sections.gjoneska_pfenning \n'
                             'goatools/test_data/sections/gjoneska_pfenning.py \n'
                             'data/gjoneska_pfenning/sections_in.txt\n'))
        p.add_argument('--outfile_detail', type=str,
                       help=('Write enrichment results into a text file \n'
                             'containing the following information: \n'
                             '1) GOEA GO terms, grouped into sections \n\n'
                             '2) List of genes and ASCII art showing section membership \n'
                             '3) Detailed list of each gene and GO terms w/their P-values \n'))
        p.add_argument('--compare', dest='compare', default=False,
                       action='store_true',
                       help="the population file as a comparison group. if this "
                       "flag is specified, the population is used as the study "
                       "plus the `population/comparison`")
        p.add_argument('--ratio', dest='ratio', type=float, default=None,
                       help="only show values where the difference between study "
                       "and population ratios is greater than this. useful for "
                       "excluding GO categories with small differences, but "
                       "containing large numbers of genes. should be a value "
                       "between 1 and 2. ")
        p.add_argument('--indent', dest='indent', default=False,
                       action='store_true', help="indent GO terms")
        p.add_argument('--obo', default="go-basic.obo", type=str,
                       help="Specifies location and name of the obo file")
        p.add_argument('--no_propagate_counts', default=False, action='store_true',
                       help="Do not propagate counts to parent terms")
        # no -r:   args.relationship == False
        # -r seen: args.relationship == True
        p.add_argument('-r', '--relationship', action='store_true',
                       help='Propagate counts up all relationships')
        # NO --relationships                -> None
        # --relationships part_of regulates -> relationships=['part_of', 'regulates']
        # --relationships=part_of           -> relationships=['part_of']
        # --relationships=part_of,regulates -> relationships=['part_of', 'regulates']
        # --relationships=part_of regulates -> NOT VALID
        p.add_argument('--relationships', nargs='*',
                       help='Propagate counts up user-specified relationships')
        p.add_argument('--method', default="bonferroni,sidak,holm,fdr_bh", type=str,
                       help=Methods().getmsg_valid_methods())
        p.add_argument('--pvalcalc', default="fisher", type=str,
                       help=str(FisherFactory()))
        p.add_argument('--min_overlap', default=0.7, type=float,
                       help="Check that a minimum amount of study genes are in the population")
        p.add_argument('--goslim', default='goslim_generic.obo', type=str,
                       help="The GO slim file is used when grouping GO terms.")
        p.add_argument('--ev_inc', type=str,
                       help="Include specified evidence codes and groups separated by commas")
        p.add_argument('--ev_exc', type=str,
                       help="Exclude specified evidence codes and groups separated by commas")
        p.add_argument('--ev_help', dest='ev_help', action='store_false',
                       help="Print all Evidence codes, with descriptions")
        p.add_argument('--ev_help_short', dest='ev_help_short', action='store_false',
                       help="Print all Evidence codes")
        # remove_goids: TBD
        #   None (Default) Remove a small number (~14) of broad GO IDs from the association
        #   True           Remove a slightly larger number of broad GO IDs (~100)
        #   False          Do not remove any broad GO IDs
        ## p.add_argument('--remove_goids', dest='remove_goids', default=None,
        ##                help="User-specified list of broad GO IDs to remove")

        if len(sys.argv) == 1:
            sys.exit(not p.print_help())
        self._prt_evidence_codes(set(sys.argv[1:]))
        args = p.parse_args()  # Namespace object from argparse
        self._adjust_relationships(args)
        self._check_input_files(args, p)
        return args

    @staticmethod
    def _prt_evidence_codes(args):
        if not {'--ev_help', '--ev_help_short'}.isdisjoint(args):
            print('\nEVIDENCE CODE HELP: --ev_exc --ev_inc')
            print('Use any of these group names, ')
            print('like Experimental or Similarity or Experimental,Similarity,')
            print('or evidence codes, like IEA or ISS,ISO,ISA in --ev_exc or --ev_inc:')
            obj = EvidenceCodes()
            if '--ev_help' in args:
                print('')
                obj.prt_details()
            if '--ev_help_short' in args:
                print('')
                obj.prt_summary_code()
            sys.exit(0)

    @staticmethod
    def _check_input_files(nspc, parser):
        """check filename args. otherwise if one of the 3 filenames is bad
        it's hard to tell which one"""
        if not len(nspc.filenames) == 3:
            parser.print_help()
            msg = """
      3 Expected files; Expected content: study population association",
      {} Actual   files: {}""".format(len(nspc.filenames), ' '.join(nspc.filenames))
            raise Exception(msg)
        for fin in nspc.filenames:
            if not os.path.exists(fin):
                return "*{}* does not exist".format(fin)
        return False

    @staticmethod
    def _adjust_relationships(args):
        """Adjust relationships for various user input"""
        # NO --relationships                -> None
        # --relationships part_of regulates -> relationships=['part_of', 'regulates']
        # --relationships=part_of           -> relationships=['part_of']
        # --relationships=part_of,regulates -> relationships=['part_of,regulates']
        # --relationships=part_of regulates -> NOT VALID
        if args.relationship:
            args.relationships = RELATIONSHIP_SET
        if args.relationships is not None:
            if len(args.relationships) == 1 and ',' in args.relationships[0]:
                args.relationships = args.relationships[0].split(',')
            args.relationships = set(args.relationships)
            chk_relationships(args.relationships)



class GoeaCliFnc:
    """For running a GOEA on the command-line."""

    # pylint: disable=too-many-instance-attributes
    def __init__(self, args):
        self.args = args
        self.sections = read_sections(self.args.sections) if self.args.sections else None
        godag_optional_attrs = self._get_optional_attrs()
        self.godag = GODag(obo_file=self.args.obo, optional_attrs=godag_optional_attrs)
        ## print('ARGS GoeaCliFnc ', self.args)
        # GET: Gene2GoReader, GafReader, GpadReader, or IdToGosReader
        self.objanno = self._get_objanno(self.args.filenames[2])
        _study, _pop = self.rd_files(*self.args.filenames[:2])
        if not self.args.compare:
            # Compare population and study gene product sets
            self.chk_genes(_study, _pop, self.objanno.associations)
        self.methods = self.args.method.split(",")
        self.itemid2name = self._init_itemid2name()
        # Get GOEnrichmentStudyNS
        self.objgoeans = self._init_objgoeans(_pop)
        # Run GOEA
        self.results_all = self.objgoeans.run_study(_study)
        # Prepare for grouping, if user-specified. Create GroupItems
        self.prepgrp = GroupItems(self, self.godag.version) if self.sections else None

    def _get_objanno(self, assoc_fn):
        """Get an annotation object"""
        # Determine annotation file format from filename, if possible
        anno_type = get_anno_desc(assoc_fn, None)
        # Default annotation file format is id2gos
        if anno_type is None:
            anno_type = self.args.annofmt if self.args.annofmt else 'id2gos'
        # kws: namespaces taxid godag
        kws = self._get_kws_objanno(anno_type)
        return get_objanno(assoc_fn, anno_type, **kws)

    def _get_ns(self):
        """Return namespaces."""
        exp_nss = {'BP', 'MF', 'CC'}
        act_nss = set(self.args.ns.split(','))
        assert not act_nss.difference(exp_nss), 'EXPECTED NAMESPACES({E}); GOT({A})'.format(
            E=','.join(exp_nss), A=','.join(act_nss.difference(exp_nss)))
        return None if act_nss == exp_nss else act_nss

    def _get_kws_objanno(self, anno_type):
        """Get keyword-args for creating an Annotation object"""
        kws = {'namespaces': self._get_ns(), 'godag': self.godag}
        if anno_type == 'gene2go':
            kws['taxid'] = self.args.taxid
        return kws

    def _init_itemid2name(self):
        """Print gene symbols instead of gene IDs, if provided."""
        if not hasattr(self.args, 'id2sym'):
            return None
        fin_id2sym = self.args.id2sym
        if fin_id2sym is not None and os.path.exists(fin_id2sym):
            id2sym = {}
            cmpl = re.compile(r'^\s*(\S+)[\s,;]+(\S+)')
            with open(fin_id2sym) as ifstrm:
                for line in ifstrm:
                    mtch = cmpl.search(line)
                    if mtch:
                        id2sym[mtch.group(1)] = mtch.group(2)
            return id2sym
        return None

    def prt_results(self, goea_results):
        """Print GOEA results to the screen or to a file."""
        # objaart = self.prepgrp.get_objaart(goea_results) if self.prepgrp is not None else None
        if self.args.outfile is None:
            self._prt_results(goea_results)
        else:
            # Users can print to both tab-separated file and xlsx file in one run.
            outfiles = self.args.outfile.split(",")
            grpwr = self.prepgrp.get_objgrpwr(goea_results) if self.prepgrp else None
            if grpwr is None:
                self.prt_outfiles_flat(goea_results, outfiles)
            else:
                grpwr.prt_outfiles_grouped(outfiles)

    def prt_outfiles_flat(self, goea_results, outfiles):
        """Write to outfiles."""
        kws = {'indent':self.args.indent, 'itemid2name':self.itemid2name}
        for outfile in outfiles:
            if outfile.endswith(".xlsx"):
                self.objgoeans.wr_xlsx(outfile, goea_results, **kws)
            #elif outfile.endswith(".txt"):  # TBD
            #    pass
            else:
                self.objgoeans.wr_tsv(outfile, goea_results, **kws)

    def _prt_results(self, goea_results):
        """Print GOEA results to the screen."""
        min_ratio = self.args.ratio
        if min_ratio is not None:
            assert 1 <= min_ratio <= 2
        results_adj = OBJPRTRES.get_adj_records(goea_results, min_ratio, self.args.pval)
        OBJPRTRES.print_date(min_ratio=min_ratio, pval=self.args.pval)
        if results_adj:
            if not self.prepgrp:
                OBJPRTRES.print_results_adj(results_adj, indent=self.args.indent)
            else:
                grpwr = self.prepgrp.get_objgrpwr(results_adj)
                grpwr.prt_txt(sys.stdout)

    def get_results(self):
        """Given all GOEA results, return the significant results (< pval)."""
        return self.get_results_sig() if self.args.pval != -1.0 else self.results_all

    def _init_objgoeans(self, pop):
        """Run gene ontology enrichment analysis (GOEA)."""
        ns2assoc = self.objanno.get_ns2assc(**self._get_anno_kws())
        ## BROAD rm_goids = self._get_remove_goids()
        rm_goids = False  # BROAD
        return GOEnrichmentStudyNS(pop, ns2assoc, self.godag,
                                   propagate_counts=not self.args.no_propagate_counts,
                                   relationships=self.args.relationships,
                                   alpha=self.args.alpha,
                                   pvalcalc=self.args.pvalcalc,
                                   methods=self.methods,
                                   remove_goids=rm_goids)
    def _get_anno_kws(self):
        """Return keyword options to obtain id2gos"""
        kws = {}
        if self.args.ev_inc is not None:
            kws['ev_include'] = set(self.args.ev_inc.split(','))
        if self.args.ev_exc is not None:
            kws['ev_exclude'] = set(self.args.ev_exc.split(','))
        return kws

    def chk_genes(self, study, pop, ntsassoc=None):
        """Compare population and study gene product sets"""
        if len(pop) < len(study):
            exit("\nERROR: The study file contains more elements than the population file. "
                 "Please check that the study file is a subset of the population file.\n")
        # check the fraction of genomic ids that overlap between study and population
        overlap = self.get_overlap(study, pop)
        if overlap < 0.95:
            sys.stderr.write("\nWARNING: only {} fraction of genes/proteins in study are found in "
                             "the population background.\n\n".format(overlap))
        if overlap <= self.args.min_overlap:
            exit("\nERROR: only {} of genes/proteins in the study are found in the "
                 "background population. Please check.\n".format(overlap))
        # Population and associations
        if ntsassoc is not None:
            assc_ids = set(nt.DB_ID for nt in ntsassoc)
            if pop.isdisjoint(assc_ids):
                if self.objanno.name == 'gene2go':
                    err = ('**FATAL: NO POPULATION ITEMS SEEN IN THE NCBI gene2go ANNOTATIONS '
                           'FOR taxid({T}). TRY: --taxid=<taxid number>')
                    exit(err.format(T=next(iter(self.objanno.taxid2asscs.keys()))))
                else:
                    exit('**FATAL: NO POPULATION ITEMS SEEN IN THE ANNOTATIONS')

    def get_results_sig(self):
        """Get significant results."""
        # Only print results when uncorrected p-value < this value.
        print("{N:7,} of {M:,} results have uncorrected P-values <= {PVAL}=pval\n".format(
            N=sum(1 for r in self.results_all if r.p_uncorrected < self.args.pval),
            M=len(self.results_all),
            PVAL=self.args.pval))
        pval_fld = self.get_pval_field()
        results = [r for r in self.results_all if getattr(r, pval_fld) <= self.args.pval]
        return results

    @staticmethod
    def get_overlap(study, pop):
        """Get he ratio of study genes which are in the population."""
        return float(len(study & pop)) / len(study)

    def get_pval_field(self):
        """Get 'p_uncorrected' or the user-specified field for determining significant results."""
        pval_fld = self.args.pval_field
        # If --pval_field [VAL] was specified
        if pval_fld is not None:
            if pval_fld[:2] != 'p_':
                pval_fld = 'p_' + pval_fld
        # If only one method was used, use that instead of the uncorrected pvalue
        elif len(self.methods) == 1:
            pval_fld = 'p_' + self.methods[0]
        # Use 'uncorrected pvalue' if there are many methods & none chosen using --pval_field
        else:
            pval_fld = 'p_uncorrected'
        if self.results_all:
            assert hasattr(next(iter(self.results_all)), pval_fld), \
                'NO PVAL({P}). EXPECTED ONE OF: {E}'.format(
                    P=self.args.pval_field,
                    E=" ".join([k for k in dir(next(iter(self.results_all))) if k[:2] == 'p_']))
        return pval_fld

    def rd_files(self, study_fn, pop_fn):
        """Read files and return study and population."""
        study, pop = self._read_geneset(study_fn, pop_fn)
        print("Study: {0} vs. Population {1}\n".format(len(study), len(pop)))
        return study, pop

    def _read_geneset(self, study_fn, pop_fn):
        """Open files containing genes. Return study genes and population genes."""
        pop = set(_.strip() for _ in open(pop_fn) if _.strip())
        study = frozenset(_.strip() for _ in open(study_fn) if _.strip())
        if next(iter(pop)).isdigit():
            pop = set(int(g) for g in pop)
            study = frozenset(int(g) for g in study)
        # some times the pop is a second group to compare, rather than the
        # population in that case, we need to make sure the overlapping terms
        # are removed first
        if self.args.compare:
            common = pop & study
            pop |= study
            pop -= common
            study -= common
            sys.stderr.write("removed %d overlapping items\n" % (len(common)))
            sys.stderr.write("Set 1: {0}, Set 2: {1}\n".format(
                len(study), len(pop)))
        return study, pop

    def _get_optional_attrs(self):
        """Given keyword args, return optional_attributes to be loaded into the GODag."""
        if self.args.relationship:
            return {'relationship',}
        if self.args.relationships is not None:
            return {'relationship',}
        if self.sections:
            return {'relationship',}
        return None

    def _get_remove_goids(self):
        """Get arguments to get a list of broad GO IDs to remove"""
        # None: (Default) Remove a small number (~14) of broad GO IDs from the association
        if self.args.remove_goids is None:
            return self.args.remove_goids
        # True: Remove a slightly larger number of broad GO IDs (~100)
        if self.args.remove_goids.lower() == 'true':
            return True
        # False: Do not remove any broad GO IDs
        if self.args.remove_goids.lower() == 'false':
            return False
        if os.path.exists(self.args.remove_goids):
            GetGOs.rdtxt_gos(self.args.remove_goids, sys.stdout)
        if isinstance(self.args.remove_goids, str) and 'GO' in self.args.remove_goids:
            return self.args.remove_goids.split(',')
        return None

class GroupItems:
    """Prepare for grouping, if specified by the user."""

    def __init__(self, objcli, godag_version):
        # _goids = set(o.id for o in godag.values() if not o.children)
        _goids = set(r.GO for r in objcli.results_all)
        _tobj = TermCounts(objcli.godag, objcli.objgoeans.get_assoc())
        # pylint: disable=line-too-long
        self.gosubdag = GoSubDag(_goids, objcli.godag, relationships=True, tcntobj=_tobj, prt=sys.stdout)
        self.grprdflt = GrouperDflts(self.gosubdag, objcli.args.goslim)
        self.hdrobj = HdrgosSections(self.grprdflt.gosubdag, self.grprdflt.hdrgos_dflt, objcli.sections)
        self.pval_fld = objcli.get_pval_field()  # primary pvalue of interest
        self.ver_list = [godag_version,
                         self.grprdflt.ver_goslims,
                         "Sections: {S}".format(S=objcli.args.sections)]
        # self.objaartall = self._init_objaartall()

    def get_objgrpwr(self, goea_results):
        """Get a GrpWr object to write grouped GOEA results."""
        sortobj = self.get_sortobj(goea_results)
        return GrpWr(sortobj, self.pval_fld, ver_list=self.ver_list)

    def get_sortobj(self, goea_results, **kws):
        """Return a Grouper object, given a list of GOEnrichmentRecord."""
        nts_goea = MgrNtGOEAs(goea_results).get_goea_nts_prt(**kws)
        goids = set(nt.GO for nt in nts_goea)
        go2nt = {nt.GO:nt for nt in nts_goea}
        grprobj = Grouper("GOEA", goids, self.hdrobj, self.grprdflt.gosubdag, go2nt=go2nt)
        grprobj.prt_summary(sys.stdout)
        # hdrgo_prt", "section_prt", "top_n", "use_sections"
        sortobj = Sorter(grprobj, section_sortby=lambda nt: getattr(nt, self.pval_fld))
        return sortobj

    # @staticmethod
    # def get_objaart(goea_results, **kws):
    #     """Return a AArtGeneProductSetsOne object."""
    #     nts_goea = MgrNtGOEAs(goea_results).get_goea_nts_prt(**kws)
    #     # objaart = AArtGeneProductSetsOne(name, goea_nts, self)

    def _init_objaartall(self):
        """Get background database info for making ASCII art."""
        kws = {
            'sortgo':lambda nt: [nt.NS, nt.dcnt],
            # fmtgo=('{p_fdr_bh:8.2e} {GO} '
            # Formatting for GO terms in grouped GO list
            'fmtgo':('{hdr1usr01:2} {NS} {GO} {s_fdr_bh:8} '
                     '{dcnt:5} {childcnt:3} R{reldepth:02} '
                     '{D1:5} {GO_name} ({study_count} study genes)\n'),
            # Formatting for GO terms listed under each gene
            'fmtgo2':('{hdr1usr01:2} {NS} {GO} {s_fdr_bh:8} '
                      '{dcnt:5} R{reldepth:02} '
                      '{GO_name} ({study_count} study genes)\n'),
            # itemid2name=ensmusg2symbol}
            }
        return AArtGeneProductSetsAll(self.grprdflt, self.hdrobj, **kws)

class GrpWr:
    """Write GO term GOEA information, grouped."""

    objprtfmt = PrtFmt()

    def __init__(self, sortobj, pval_fld, ver_list):
        self.sortobj = sortobj
        self.pval_fld = pval_fld
        self.ver_list = ver_list
        self.flds_all = next(iter(self.sortobj.grprobj.go2nt.values()))._fields
        self.flds_cur = self._init_flds_cur()
        self.desc2nts = self.sortobj.get_desc2nts(hdrgo_prt=False)

    def prt_outfiles_grouped(self, outfiles):
        """Write to outfiles."""
        for outfile in outfiles:
            if outfile.endswith(".xlsx"):
                self.wr_xlsx(outfile)
            elif outfile.endswith(".txt"):
                self.wr_txt(outfile)
            else:
                self.wr_tsv(outfile)

    def wr_xlsx(self, fout_xlsx):
        """Print grouped GOEA results into an xlsx file."""
        objwr = WrXlsxSortedGos("GOEA", self.sortobj)
        #### fld2fmt['ratio_in_study'] = '{:>8}'
        #### fld2fmt['ratio_in_pop'] = '{:>12}'
        #### ntfld2wbfmtdict = {
        # ntfld_wbfmt = {
        #     'ratio_in_study': {'align':'right'},
        #     'ratio_in_pop':{'align':'right'}}
        kws_xlsx = {
            'title': self.ver_list,
            'fld2fmt': {f:'{:8.2e}' for f in self.flds_cur if f[:2] == 'p_'},
            #'ntfld_wbfmt': ntfld_wbfmt,
            #### 'ntval2wbfmtdict': ntval2wbfmtdict,
            #'hdrs': [],
            'prt_flds': self.flds_cur}
        objwr.wr_xlsx_nts(fout_xlsx, self.desc2nts, **kws_xlsx)

    def wr_tsv(self, fout_tsv):
        """Print grouped GOEA results into a tab-separated file."""
        with open(fout_tsv, 'w') as prt:
            kws_tsv = {
                'fld2fmt': {f:'{:8.2e}' for f in self.flds_cur if f[:2] == 'p_'},
                'prt_flds':self.flds_cur}
            prt_tsv_sections(prt, self.desc2nts['sections'], **kws_tsv)
            print("  WROTE: {TSV}".format(TSV=fout_tsv))

    def wr_txt(self, fout_txt):
        """Write to a file GOEA results in an ASCII text format."""
        with open(fout_txt, 'w') as prt:
            for line in self.ver_list:
                prt.write("{LINE}\n".format(LINE=line))
            self.prt_txt(prt)
            print("  WROTE: {TXT}".format(TXT=fout_txt))

    def prt_tsv(self, prt=sys.stdout):
        """Print an ASCII text format."""
        prtfmt = self.objprtfmt.get_prtfmt_str(self.flds_cur)
        prt.write("{FLDS}\n".format(FLDS=" ".join(self.flds_cur)))
        WrSectionsTxt.prt_sections(prt, self.desc2nts['sections'], prtfmt, secspc=True)

    def prt_txt(self, prt=sys.stdout):
        """Print an ASCII text format."""
        prtfmt = self.objprtfmt.get_prtfmt_str(self.flds_cur)
        prt.write("{FLDS}\n".format(FLDS=" ".join(self.flds_cur)))
        WrSectionsTxt.prt_sections(prt, self.desc2nts['sections'], prtfmt, secspc=True)

    def _init_flds_cur(self):
        """Choose fields to print from a multitude of available fields."""
        flds = []
        # ('GO', 'NS', 'enrichment', 'name', 'ratio_in_study', 'ratio_in_pop', 'depth',
        # 'p_uncorrected', 'p_bonferroni', 'p_sidak', 'p_holm', 'p_fdr_bh',
        # 'pop_n', 'pop_count', 'pop_items'
        # 'study_n', 'study_count', 'study_items',
        # 'is_ratio_different', 'level', 'is_obsolete',
        # 'namespace', 'reldepth', 'alt_ids', 'format_txt', 'hdr_idx',
        # 'is_hdrgo', 'is_usrgo', 'num_usrgos', 'hdr1usr01', 'alt', 'GO_name',
        # 'dcnt', 'D1', 'tcnt', 'tfreq', 'tinfo', 'childcnt', 'REL',
        # 'REL_short', 'rel', 'id')
        flds0 = ['GO', 'NS', 'enrichment', self.pval_fld, 'dcnt', 'tinfo', 'depth',
                 'ratio_in_study', 'ratio_in_pop', 'name']
        flds_p = [f for f in self.flds_all if f[:2] == 'p_' and f != self.pval_fld]
        flds.extend(flds0)
        if flds_p:
            flds.extend(flds_p)
        flds.append('study_count')
        flds.append('study_items')
        return flds


# Copyright (C) 2010-2019, H Tang et al. All rights reserved.
