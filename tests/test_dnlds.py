"""Test downloading go-basic.obo file."""

import os
from goatools.base import download_go_basic_obo
from goatools.base import download_ncbi_associations

def test_obo():
    """Test downloading of Ontology file."""
    fdnld = download_go_basic_obo()
    os.system("rm -f {FILE}".format(FILE=fdnld))
    fdnld = download_go_basic_obo()
    assert os.path.isfile(fdnld)

def test_NCBI_assc():
    """Test downloading of associations from NCBI."""
    fdnld = download_ncbi_associations()
    os.system("rm -f {FILE}".format(FILE=fdnld))
    fdnld = download_ncbi_associations()
    assert os.path.isfile(fdnld)

if __name__ == '__main__':
    test_obo()
    test_NCBI_assc()
