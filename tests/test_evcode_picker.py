#!/usr/bin/env python
"""Test that expected 3 letter codes are chosen when given: inc, exc for codes and groups"""

__copyright__ = "Copyright (C) 2016-2019, DV Klopfenstein, H Tang. All rights reserved."
__author__ = "DV Klopfenstein"

import os
from goatools.evidence_codes import EvidenceCodes

REPO = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..")


def test_evcode_picker():
    """Test that expected 3 letter codes are chosen when given: inc, exc for codes and groups"""
    obj = EvidenceCodes()
    # pylint: disable=superfluous-parens
    act = obj.get_evcodes()
    print('ALL POSITIVE CODES: {C}'.format(C=' '.join(sorted(act))))
    assert 'ND' not in act and len(act) > 15, act
    #
    act = obj.get_evcodes({'Experimental'})
    assert act == set(['EXP', 'IDA', 'IPI', 'IMP', 'IGI', 'IEP']), act
    #
    act = obj.get_evcodes({'Experimental'}, {'IEP'})
    assert act == set(['EXP', 'IDA', 'IPI', 'IMP', 'IGI']), act
    #
    act = obj.get_evcodes({'Experimental', 'Similarity'}, {'IEP', 'IMR'})
    exp = {
        'EXP', 'IDA', 'IPI', 'IMP', 'IGI',
        'ISS', 'ISO', 'ISA', 'ISM', 'IGC', 'IBA', 'IBD', 'IKR', 'IRD'}
    assert act == exp, act
    #
    act = obj.get_evcodes(None, {'IEA'})
    exp = set(obj.code2nt)
    exp.difference_update({'IEA', 'ND'})
    assert act == exp, act.symmetric_difference(exp)
    #
    obj.prt_details()
    obj.prt_summary_code()
    print("**TEST PASSED")


if __name__ == '__main__':
    test_evcode_picker()
