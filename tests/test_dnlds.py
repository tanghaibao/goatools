"""tests/test_dnlds.py: Test downloading go-basic.obo file and goslims obo/owl."""

import os
import sys
from goatools.base import download_go_basic_obo
from goatools.base import download_ncbi_associations

REPO = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..")


def test_dnlds(prt=sys.stdout):
    """Test downloads of ontologies and NCBI associations."""
    goslims = [
        'goslim_aspergillus',
        'goslim_candida',
        #'goslim_chembl',  # 404 not found
        'goslim_generic',
        'goslim_metagenomics',
        'goslim_pir',
        'goslim_plant',
        'goslim_pombe',
        'goslim_synapse',
        'goslim_virus',
        'goslim_yeast']
    # Test downloads of ontologies.
    dnld_ontology(os.path.join(REPO, "go-basic.obo"))
    # Test downloads of go-slim ontologies.
    for goslim in goslims:
        for ext in ['obo', 'owl']:
            file_dst = os.path.join(REPO, "{DAG}.{EXT}".format(DAG=goslim, EXT=ext))
            dnld_ontology(file_dst)
    # Test downloading of associations from NCBI.
    file_assc = os.path.join(REPO, "gene2go")
    #os.system("rm -f {FILE}".format(FILE=file_assc))
    #download_ncbi_associations(file_assc, prt, loading_bar=None)
    #download_ncbi_associations(file_assc, prt, loading_bar=None)
    #assert os.path.isfile(file_assc), "FILE({F}) EXPECTED TO EXIST".format(F=file_assc)

def dnld_ontology(filename):
    """Test downloading of ontologies."""
    # download_go_basic_obo(filename, loading_bar=None)
    os.system("rm -f {FILE}".format(FILE=filename))
    download_go_basic_obo(filename, loading_bar=None)
    download_go_basic_obo(filename, loading_bar=None)
    assert os.path.isfile(filename), "FILE({F}) EXPECTED TO EXIST".format(F=filename)

if __name__ == '__main__':
    test_dnlds()
