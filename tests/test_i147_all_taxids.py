#!/usr/bin/env python3
"""Work with all taxids using Gene2GoReader"""
# https://github.com/tanghaibao/goatools/issues/147

from os import system
from os.path import join
from os.path import exists
from goatools.base import download_go_basic_obo
from goatools.base import download_ncbi_associations
from goatools.obo_parser import GODag
from goatools.anno.genetogo_reader import Gene2GoReader
from tests.utils import REPO

def test_i147_all_taxids():
    """Work with all taxids using Gene2GoReader"""
    # 1. Download Ontologies and Associations
    # 1a. Download Ontologies, if necessary
    #     Get http://geneontology.org/ontology/go-basic.obo
    fin_dag = join(REPO, "go-basic.obo")
    download_go_basic_obo(fin_dag)

    # 1b. Download Associations, if necessary
    #     Get ftp://ftp.ncbi.nlm.nih.gov/gene/DATA/gene2go.gz
    fin_anno = join(REPO, "gene2go")
    if exists(fin_anno):
        system('rm -f gene2go*')
    fin_gene2go = download_ncbi_associations(fin_anno)
    assert fin_anno == fin_gene2go

    # 2. Load Ontologies, Associations and Background gene set
    # 2a. Load Ontologies
    godag = GODag(fin_dag)

    # 2b. Load Associations for all species
    #     Read NCBI's gene2go. Store annotations in a list of namedtuples
    objanno_all = Gene2GoReader(fin_gene2go, godag=godag, taxids=True)
    objanno_mmu = Gene2GoReader(fin_gene2go, godag=godag, taxids=[10090])
    objanno_mmuhsa = Gene2GoReader(fin_gene2go, godag=godag, taxids=[10090, 9606])

    # Get associations
    # pylint: disable=bad-whitespace
    ns2assoc_all_mmu = _run_get_ns2assc(10090, objanno_all)
    ns2assoc_mmu_mmu = _run_get_ns2assc(10090, objanno_mmu)
    ns2assoc_mmuhsa_all = _run_get_ns2assc(True,  objanno_mmuhsa)
    ns2assoc_mmuhsa_mmu = _run_get_ns2assc(10090, objanno_mmuhsa)

    # Check results
    for nspc in ['BP', 'MF', 'CC']:
        assert ns2assoc_mmu_mmu[nspc] == ns2assoc_all_mmu[nspc]
        assert ns2assoc_mmu_mmu[nspc] == ns2assoc_mmuhsa_mmu[nspc]
    _chk_mmuhsa_all(objanno_mmuhsa, objanno_all, ns2assoc_mmuhsa_all)

def _chk_mmuhsa_all(objanno_mmuhsa, objanno_all, ns2assoc_mmuhsa_all):
    """Check combining multiple species"""
    #   1. taxid=None returns all loaded annotations
    #   2. taxid=[10090, 9606] returns only those annotations, not all annotations
    ns2gene2gos_2of2 = objanno_mmuhsa.get_ns2assc(True)
    ns2gene2gos_2ofall = objanno_all.get_ns2assc([10090, 9606])
    _prt_assc_counts(ns2gene2gos_2of2)
    _prt_assc_counts(ns2gene2gos_2ofall)
    assert ns2gene2gos_2of2 == ns2gene2gos_2ofall
    assert ns2gene2gos_2of2 == objanno_all.get_ns2assc({10090, 9606})
    assert ns2gene2gos_2of2 == ns2assoc_mmuhsa_all
    assert objanno_all.get_ns2assc() == {}
    for nspc in ['BP', 'MF', 'CC']:
        assert ns2gene2gos_2of2[nspc] == ns2gene2gos_2ofall[nspc]

def _prt_assc_counts(ns2assc):
    """Print the number of genes and GO IDs in an association"""
    for nspc, gene2goids in sorted(ns2assc.items()):
        print("{NS} {N:6,} genes, {GOs:6,} GOs".format(
            NS=nspc, N=len(gene2goids), GOs=len(set.union(*gene2goids.values()))))

def _run_get_ns2assc(taxid, objanno):
    # Get namespace2association where:
    #    namespace is:
    #        BP: biological_process
    #        MF: molecular_function
    #        CC: cellular_component
    #    assocation is a dict:
    #        key: NCBI GeneID
    #        value: A set of GO IDs associated with that gene
    ns2assoc = objanno.get_ns2assc(taxid)

    num_nts = sum(len(a) for a in objanno.taxid2asscs.values())
    print('HEADER {T} anno taxids, {N:,} nts, {U} user-requested taxid:'.format(
        T=len(objanno.taxid2asscs), N=num_nts, U=taxid))
    for nspc, id2gos in ns2assoc.items():
        print("  DATA {NS} {N:,} annotated genes for taxids: {T}".format(
            NS=nspc, N=len(id2gos),
            T=sorted(objanno.taxid2asscs.keys()) if taxid is None else taxid))
    return ns2assoc

if __name__ == '__main__':
    test_i147_all_taxids()
