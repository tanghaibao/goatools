"""Test downloading go-basic.obo file."""

import os
from goatools.base import download_go_basic_obo
from goatools.base import download_ncbi_associations

def test_dnlds():
    """Test downloads of ontologies and NCBI associations."""
    # Test downloads of ontologies.
    cwd = os.getcwd()
    file_obo = os.path.join(cwd, "go-basic.obo")
    download_go_basic_obo(file_obo, loading_bar=None)
    os.system("rm -f {FILE}".format(FILE=file_obo))
    download_go_basic_obo(file_obo, loading_bar=None)
    assert os.path.isfile(file_obo)
    # Test downloading of associations from NCBI.
    file_assc = os.path.join(cwd, "gene2go")
    download_ncbi_associations(file_assc, loading_bar=None)
    os.system("rm -f {FILE}".format(FILE=file_assc))
    download_ncbi_associations(file_assc, loading_bar=None)
    assert os.path.isfile(file_assc)

if __name__ == '__main__':
    test_dnlds()
