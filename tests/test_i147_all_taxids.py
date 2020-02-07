#!/usr/bin/env python3
"""Work with all taxids using Gene2GoReader"""
# https://github.com/tanghaibao/goatools/issues/147

from goatools.base import download_go_basic_obo
from goatools.base import download_ncbi_associations
from goatools.obo_parser import GODag
from goatools.anno.genetogo_reader import Gene2GoReader

def test_i147_all_taxids():
    """Work with all taxids using Gene2GoReader"""
    # 1. Download Ontologies and Associations
    # 1a. Download Ontologies, if necessary
    #     Get http://geneontology.org/ontology/go-basic.obo
    download_go_basic_obo()

    # 1b. Download Associations, if necessary
    #     Get ftp://ftp.ncbi.nlm.nih.gov/gene/DATA/gene2go.gz
    fin_gene2go = download_ncbi_associations()

    # 2. Load Ontologies, Associations and Background gene set
    # 2a. Load Ontologies
    godag = GODag("go-basic.obo")

    # 2b. Load Associations for all species
    #     Read NCBI's gene2go. Store annotations in a list of namedtuples
    objanno = Gene2GoReader(fin_gene2go, godag=godag, taxids=True)

    # Get namespace2association where:
    #    namespace is:
    #        BP: biological_process
    #        MF: molecular_function
    #        CC: cellular_component
    #    assocation is a dict:
    #        key: NCBI GeneID
    #        value: A set of GO IDs associated with that gene
    ns2assoc = objanno.get_ns2assc()

    for nspc, id2gos in ns2assoc.items():
        print("{NS} {N:,} annotated mouse genes".format(NS=nspc, N=len(id2gos)))

if __name__ == '__main__':
    test_i147_all_taxids()
