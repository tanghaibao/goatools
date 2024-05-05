"""Run GOATOOLS Gene Ontology Analysis on Nature 3102 data. Return results."""

__copyright__ = "Copyright (C) 2016-2018, DV Klopfenstein, H Tang. All rights reserved."
__author__ = "DV Klopfenstein"

import os
import pandas as pd

from tests.utils import repofn

from ..anno.genetogo_reader import Gene2GoReader
from ..associations import dnld_ncbi_gene_file
from ..base import get_godag
from ..go_enrichment import GOEnrichmentStudy

from .genes_NCBI_10090_ProteinCoding import GENEID2NT as GeneID2nt_mus


def get_goea_results(keep_if=None):
    """Demonstrate printing a subset of all available fields using two methods."""
    if keep_if is None:
        keep_if = (
            lambda nt: getattr(nt, "p_fdr_bh") < 0.05
        )  # keep if results are significant
    # 1. Gene Ontology Enrichment Analysis
    #    1a. Initialize: Load ontologies, associations, and population gene IDs
    taxid = 10090  # Mouse study
    geneids_pop = GeneID2nt_mus.keys()  # Mouse protein-coding genes
    goeaobj = get_goeaobj("fdr_bh", geneids_pop, taxid)
    #    1b. Run GOEA
    geneids_study = get_geneid2symbol("nbt.3102-S4_GeneIDs.xlsx")
    return {
        "goea_results": goeaobj.run_study(geneids_study, keep_if=keep_if),
        "goeaobj": goeaobj,
        "geneids_study": geneids_study,
        "geneids_pop": geneids_pop,
        "obo_dag": goeaobj.obo_dag,
    }


def get_geneid2symbol(fin_xlsx):
    """Read xlsx file return dictionary with Entrez GeneID keys to Symbol data."""
    gene2symbol = {}
    data_dir = os.path.join(os.path.dirname(os.path.abspath(__file__)), "nbt_3102")
    tbl_genes = f"{data_dir}/{fin_xlsx}"
    df = pd.read_excel(
        tbl_genes, engine="openpyxl", header=None, names=["Symbol", "GeneID", "Desc"]
    )
    for _, row in df.iterrows():
        symbol = row["Symbol"]
        geneid = int(row["GeneID"])
        if geneid:
            gene2symbol[geneid] = symbol
    return gene2symbol


def get_goeaobj(method, geneids_pop, taxid, nspc="BP"):
    """Load: ontologies, associations, and population geneids."""
    fin_obo = os.path.join(os.getcwd(), "go-basic.obo")
    godag = get_godag(fin_obo)
    assoc_geneid2gos = get_annotations(taxid, nspc)
    goeaobj = GOEnrichmentStudy(
        geneids_pop,
        assoc_geneid2gos,
        godag,
        propagate_counts=False,
        alpha=0.05,
        methods=[method],
    )
    # godag is also found in goeaobj.godag
    return goeaobj


def get_annotations(taxid, nspc="BP"):
    """Download annotations. Return dict of gene-to-GOs"""
    fin_anno = repofn("gene2go")
    dnld_ncbi_gene_file(fin_anno)
    objanno = Gene2GoReader(fin_anno, taxid=taxid)
    return objanno.get_id2gos(namespace=nspc)


# Copyright (C) 2016-2018, DV Klopfenstein, H Tang. All rights reserved.
