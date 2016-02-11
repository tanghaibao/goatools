import sys
import os
from goatools.go_enrichment import GOEnrichmentStudy
from goatools.obo_parser import GODag
from goatools.associations import read_associations

"""Test that GOEnrichmentStudy fails elegantly given incorrect stimulus.

        python test_goea_errors.py
"""

__copyright__ = "Copyright (C) 2016, DV Klopfenstein, H Tang. All rights reserved."
__author__ = "DV Klopfenstein"

ROOT = os.path.dirname(os.path.abspath(__file__)) + "/../data/"


def init_goea(**kws):
    """Initialize GODag and GOEnrichmentStudy."""
    obo_dag = GODag(ROOT + "go-basic.obo")
    assoc = read_associations(ROOT + "association", no_top=True)
    popul_ids = [line.rstrip() for line in open(ROOT + "population")]
    methods = kws['methods'] if 'methods' in kws else ['not_bonferroni']
    study_ids = [line.rstrip() for line in open(ROOT + "study")]
    return GOEnrichmentStudy(popul_ids, assoc, obo_dag, methods=methods), study_ids


def run_method_bad_ini():
    """Test attempting to use an unsupported method in initialization."""
    goea, study_ids = init_goea(methods=['not_fdr'])
    # Test that method(s) set during initialization are valid
    goea.run_study(study_ids)


def run_method_bad_run():
    """Test attempting to use an unsupported method in run."""
    goea, study_ids = init_goea()
    # Test that method(s) set while running a GOEA on a study are valid
    goea.run_study(study_ids, methods=['invalid_method'])


def test_all(log=sys.stdout):
    """Run all tests."""
    tests = [
        (run_method_bad_ini, "INVALID METHOD(not_fdr)"),
        (run_method_bad_run, "INVALID METHOD(invalid_method)"),
    ]
    for test, exp_errmsg in tests:
        try:
            test()
        except Exception as inst:
            # Run next test
            if str(inst).startswith(exp_errmsg):
                log.write("Test PASSED. Expected error message seen: {EXP}\n".format(
                    EXP=exp_errmsg))
            else:
                raise Exception("EXPECTED({EXP}). ACTUAL({ACT})".format(
                    EXP=exp_errmsg, ACT=inst))


if __name__ == '__main__':
    test_all()

# Copyright (C) 2016, DV Klopfenstein, H Tang. All rights reserved.
