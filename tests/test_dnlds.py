"""tests/test_dnlds.py: Test downloading go-basic.obo file and goslims obo/owl."""

import os

from goatools.base import download_go_basic_obo

REPO = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..")


def test_dnlds():
    """Test downloads of ontologies and NCBI associations."""
    goslims = [
        "goslim_aspergillus",
        "goslim_candida",
        #'goslim_chembl',  # 404 not found
        "goslim_generic",
        "goslim_metagenomics",
        "goslim_pir",
        "goslim_plant",
        "goslim_pombe",
        "goslim_synapse",
        "goslim_virus",
        "goslim_yeast",
    ]
    # Test downloads of ontologies.
    dnld_ontology(os.path.join(REPO, "go-basic.obo"))
    # Test downloads of go-slim ontologies.
    for goslim in goslims:
        for ext in ["obo", "owl"]:
            file_dst = os.path.join(REPO, "{DAG}.{EXT}".format(DAG=goslim, EXT=ext))
            dnld_ontology(file_dst)


def dnld_ontology(filename):
    """Test downloading of ontologies."""
    os.system("rm -f {FILE}".format(FILE=filename))
    download_go_basic_obo(filename)
    download_go_basic_obo(filename)
    assert os.path.isfile(filename), "FILE({F}) EXPECTED TO EXIST".format(F=filename)


if __name__ == "__main__":
    test_dnlds()
