"""
Routines to read in association file between genes and GO terms.
"""

__copyright__ = "Copyright (C) 2010-2019, H Tang et al. All rights reserved."
__author__ = "various"

from collections import defaultdict
import os
import sys
from goatools.base import dnld_file
from goatools.base import ftp_get
from goatools.anno.factory import get_objanno
from goatools.anno.factory import get_anno_desc
from goatools.anno.factory import get_objanno_g_kws
from goatools.semantic import TermCounts
from goatools.anno.gaf_reader import GafReader
from goatools.anno.genetogo_reader import Gene2GoReader
from goatools.anno.opts import AnnoOptions
from goatools.utils import get_b2aset as utils_get_b2aset

def dnld_assc(assc_name, go2obj=None, namespace='BP', prt=sys.stdout):
    """Download association from http://geneontology.org/gene-associations."""
    # Example assc_name: "tair.gaf"
    # Download the Association
    dirloc, assc_base = os.path.split(assc_name)
    if not dirloc:
        dirloc = os.getcwd()
    assc_locfile = os.path.join(dirloc, assc_base) if not dirloc else assc_name
    dnld_annotation(assc_locfile, prt)
    # Read the downloaded nt120GV)association
    assc_orig = read_gaf(assc_locfile, namespace=namespace, godag=go2obj, prt=prt)
    if go2obj is None:
        return assc_orig
    # If a GO DAG is provided, use only GO IDs present in the GO DAG
    assc = {}
    goids_dag = set(go2obj.keys())
    for gene, goids_cur in assc_orig.items():
        assc[gene] = goids_cur.intersection(goids_dag)
    return assc

def dnld_annotation(assc_file, prt=sys.stdout):
    """Download gaf, gpad, or gpi from http://current.geneontology.org/annotations/"""
    if not os.path.isfile(assc_file):
        # assc_http = "http://geneontology.org/gene-associations/"
        assc_http = "http://current.geneontology.org/annotations/"
        _, assc_base = os.path.split(assc_file)
        src = os.path.join(assc_http, "{ASSC}.gz".format(ASSC=assc_base))
        dnld_file(src, assc_file, prt, loading_bar=None)

def read_associations(assoc_fn, anno_type='id2gos', namespace='BP', **kws):
    """Return associatinos in id2gos format"""
    # kws get_objanno: taxids hdr_only prt allow_missing_symbol
    obj = get_objanno(assoc_fn, anno_type, **kws)
    # kws get_id2gos: ev_include ev_exclude keep_ND keep_NOT b_geneid2gos go2geneids
    return obj.get_id2gos(namespace, **kws)

def get_assoc_ncbi_taxids(taxids, force_dnld=False, loading_bar=True, **kws):
    """Download NCBI's gene2go. Return annotations for user-specified taxid(s)."""
    print('DEPRECATED read_ncbi_gene2go: USE Gene2GoReader FROM goatools.anno.genetogo_reader')
    # pylint: disable=protected-access
    frm = sys._getframe().f_back.f_code
    print('DEPRECATED read_ncbi_gene2go CALLED FROM: {PY} BY {FNC}'.format(
        PY=frm.co_filename, FNC=frm.co_name))
    fin = kws['gene2go'] if 'gene2go' in kws else os.path.join(os.getcwd(), "gene2go")
    dnld_ncbi_gene_file(fin, force_dnld, loading_bar=loading_bar)
    return read_ncbi_gene2go(fin, taxids, **kws)

# pylint: disable=unused-argument
def dnld_ncbi_gene_file(fin, force_dnld=False, log=sys.stdout, loading_bar=True):
    """Download a file from NCBI Gene's ftp server."""
    if not os.path.exists(fin) or force_dnld:
        import gzip
        fin_dir, fin_base = os.path.split(fin)
        fin_gz = "{F}.gz".format(F=fin_base)
        fin_gz = os.path.join(fin_dir, fin_gz)
        if os.path.exists(fin_gz):
            os.remove(fin_gz)
        fin_ftp = "ftp://ftp.ncbi.nlm.nih.gov/gene/DATA/{F}.gz".format(F=fin_base)
        ## if log is not None:
        ##     log.write("  DOWNLOADING GZIP: {GZ}\n".format(GZ=fin_ftp))
        ## if loading_bar:
        ##     loading_bar = wget.bar_adaptive
        ## wget.download(fin_ftp, bar=loading_bar)
        ## rsp = wget(fin_ftp)
        ftp_get(fin_ftp, fin_gz)
        with gzip.open(fin_gz, 'rb') as zstrm:
            if log is not None:
                log.write("\n  READ GZIP:  {F}\n".format(F=fin_gz))
            with open(fin, 'wb') as ostrm:
                ostrm.write(zstrm.read())
                if log is not None:
                    log.write("  WROTE UNZIPPED: {F}\n".format(F=fin))

def dnld_annofile(fin_anno, anno_type):
    """Download annotation file, if needed"""
    if os.path.exists(fin_anno):
        return
    anno_type = get_anno_desc(fin_anno, anno_type)
    if anno_type == 'gene2go':
        dnld_ncbi_gene_file(fin_anno)
    if anno_type in {'gaf', 'gpad'}:
        dnld_annotation(fin_anno)

def read_ncbi_gene2go(fin_gene2go, taxids=None, namespace='BP', **kws):
    """Read NCBI's gene2go. Return gene2go data for user-specified taxids."""
    print('DEPRECATED read_ncbi_gene2go: USE Gene2GoReader FROM goatools.anno.genetogo_reader')
    # pylint: disable=protected-access
    frm = sys._getframe().f_back.f_code
    print('DEPRECATED read_ncbi_gene2go CALLED FROM: {PY} BY {FNC}'.format(
        PY=frm.co_filename, FNC=frm.co_name))
    obj = Gene2GoReader(fin_gene2go, taxids=taxids)
    # By default, return id2gos. User can cause go2geneids to be returned by:
    #   >>> read_ncbi_gene2go(..., go2geneids=True
    if 'taxid2asscs' not in kws:
        if len(obj.taxid2asscs) == 1:
            taxid = next(iter(obj.taxid2asscs))
            kws_ncbi = {k:v for k, v in kws.items() if k in AnnoOptions.keys_exp}
            kws_ncbi['taxid'] = taxid
            return obj.get_id2gos(namespace, **kws_ncbi)
    # Optional detailed associations split by taxid and having both ID2GOs & GO2IDs
    # e.g., taxid2asscs = defaultdict(lambda: defaultdict(lambda: defaultdict(set))
    t2asscs_ret = obj.get_taxid2asscs(taxids, **kws)
    t2asscs_usr = kws.get('taxid2asscs', defaultdict(lambda: defaultdict(lambda: defaultdict(set))))
    if 'taxid2asscs' in kws:
        obj.fill_taxid2asscs(t2asscs_usr, t2asscs_ret)
    return obj.get_id2gos_all(t2asscs_ret)

def get_gaf_hdr(fin_gaf):
    """Read Gene Association File (GAF). Return GAF version and data info."""
    return GafReader(fin_gaf, hdr_only=True).hdr

# pylint: disable=line-too-long
def read_gaf(fin_gaf, prt=sys.stdout, hdr_only=False, namespace='BP', allow_missing_symbol=False, **kws):
    """Read Gene Association File (GAF). Return data."""
    return GafReader(
        fin_gaf, hdr_only=hdr_only, prt=prt, allow_missing_symbol=allow_missing_symbol, godag=kws.get('godag')).get_id2gos(
            namespace, **kws)

def get_b2aset(a2bset):
    """Given gene2gos, return go2genes. Given go2genes, return gene2gos."""
    print('DEPRECATED get_b2aset MOVED: USE get_b2aset IN goatools.utils')
    # pylint: disable=protected-access
    frm = sys._getframe().f_back.f_code
    print('DEPRECATED get_b2aset CALLED FROM: {PY} BY {FNC}'.format(PY=frm.co_filename, FNC=frm.co_name))
    return utils_get_b2aset(a2bset)

def get_assc_pruned(assc_geneid2gos, min_genecnt=None, max_genecnt=None, prt=sys.stdout):
    """Remove GO IDs associated with large numbers of genes. Used in stochastic simulations."""
    # DEFN WAS: get_assc_pruned(assc_geneid2gos, max_genecnt=None, prt=sys.stdout):
    #      ADDED min_genecnt argument and functionality
    if max_genecnt is None and min_genecnt is None:
        return assc_geneid2gos, set()
    go2genes_orig = utils_get_b2aset(assc_geneid2gos)
    # go2genes_prun = {go:gs for go, gs in go2genes_orig.items() if len(gs) <= max_genecnt}
    go2genes_prun = {}
    for goid, genes in go2genes_orig.items():
        num_genes = len(genes)
        if (min_genecnt is None or num_genes >= min_genecnt) and \
           (max_genecnt is None or num_genes <= max_genecnt):
            go2genes_prun[goid] = genes
    num_was = len(go2genes_orig)
    num_now = len(go2genes_prun)
    gos_rm = set(go2genes_orig.keys()).difference(set(go2genes_prun.keys()))
    assert num_was-num_now == len(gos_rm)
    if prt is not None:
        if min_genecnt is None:
            min_genecnt = 1
        if max_genecnt is None:
            max_genecnt = "Max"
        prt.write("{N:4} GO IDs pruned. Kept {NOW} GOs assc w/({m} to {M} genes)\n".format(
            m=min_genecnt, M=max_genecnt, N=num_was-num_now, NOW=num_now))
    return utils_get_b2aset(go2genes_prun), gos_rm

def read_annotations(**kws):
    """Read annotations from either a GAF file or NCBI's gene2go file."""
    # Read and save annotation lines
    objanno = get_objanno_g_kws(**kws)
    # Return associations
    return objanno.get_id2gos() if objanno is not None else {}

def get_tcntobj(go2obj, **kws):
    """Return a TermCounts object if the user provides an annotation file, otherwise None."""
    # kws: gpad gaf gene2go id2gos
    objanno = get_objanno_g_kws(**kws)
    if objanno:
        return TermCounts(go2obj, objanno.get_id2gos_nss())
    return None


# Copyright (C) 2010-2019, H Tang et al. All rights reserved."
