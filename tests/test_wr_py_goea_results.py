"""_resultsTest writing GOATOOLS Gene Ontology Gene Enrichment results to a Python module."""

import os
import importlib
from goatools.test_data.nature3102_goea import get_goea_results
from goatools.nt_utils import wr_py_nts
from goatools.rpt.goea_nt_xfrm import get_goea_nts_prt

def test_wrpy():
    """Test writing GOATOOLS GOEA results to a Python module as a list of nts."""
    # 1. Run GOATOOLS Gene Ontology Enrichment Analysis
    nature_data = get_goea_results()
    # 2. Convert GOATOOLS GOEA results into a list of namedtuples
    goea_results = nature_data['goea_results']
    nts_goea = get_goea_nts_prt(goea_results)

    # 3. Save GOATOOLS GOEA into a Python module
    #    3a. Python module name
    module_name = "nbt3102_goea"
    fout_py = get_fout_py(module_name)
    #    3b. Save GOATOOLS GOEA into a Python module
    wr_py_nts(fout_py, nts_goea, varname="nts")
    # 4. Check results
    nts_py = importlib.import_module(module_name).nts
    assert len(nts_goea) == len(nts_py)

    # Alternatively, save module through goea object
    module_name = "nbt3102_goea_alt"
    fout_py = get_fout_py(module_name)
    nature_data['goeaobj'].wr_py_goea_results(fout_py, goea_results)
    nts_py = importlib.import_module(module_name).goea_results
    assert len(nts_goea) == len(nts_py)

def get_fout_py(module_name):
    """Given a module name, return the name of the corresponding Python file."""
    fout_py = "{MODULE}.py".format(MODULE=module_name)
    if os.path.isfile(fout_py):
        os.remove(fout_py)
    return fout_py


if __name__ == '__main__':
    test_wrpy()
