"""Test an object of class, Methods."""

import collections as cx
from goatools.multiple_testing import Methods

def test_init_methods():
    """Test initializing methods."""
    mobj = Methods()
    assert mobj._srcmethod2fieldname == get_exp_fieldnames()
    assert mobj.getmsg_valid_methods() == get_expstr_fieldnames()
    assert mobj.methods == [mobj.NtMethodInfo(source='local', method='bonferroni', fieldname='bonferroni')]
    mobj._add_method_src('statsmodels', 'fdr_bh')
    assert mobj.methods == [
        mobj.NtMethodInfo(source='local', method='bonferroni', fieldname='bonferroni'), 
        mobj.NtMethodInfo(source='statsmodels', method='fdr_bh', fieldname='fdr_bh')]
    sm_methods = ['sm_{}'.format(m) for m in mobj.all_methods[1][1]] # statsmodels
    mobj._init_methods(sm_methods)
    assert mobj.methods == [
        mobj.NtMethodInfo(source='statsmodels', method='bonferroni', fieldname='sm_bonferroni'), 
        mobj.NtMethodInfo(source='statsmodels', method='sidak', fieldname='sm_sidak'), 
        mobj.NtMethodInfo(source='statsmodels', method='holm-sidak', fieldname='sm_holm-sidak'), 
        mobj.NtMethodInfo(source='statsmodels', method='holm', fieldname='sm_holm'), 
        mobj.NtMethodInfo(source='statsmodels', method='simes-hochberg', fieldname='sm_simes-hochberg'), 
        mobj.NtMethodInfo(source='statsmodels', method='hommel', fieldname='sm_hommel'), 
        mobj.NtMethodInfo(source='statsmodels', method='fdr_bh', fieldname='sm_fdr_bh'), 
        mobj.NtMethodInfo(source='statsmodels', method='fdr_by', fieldname='sm_fdr_by'), 
        mobj.NtMethodInfo(source='statsmodels', method='fdr_tsbh', fieldname='sm_fdr_tsbh'), 
        mobj.NtMethodInfo(source='statsmodels', method='fdr_tsbky', fieldname='sm_fdr_tsbky'), 
        mobj.NtMethodInfo(source='statsmodels', method='fdr_gbs', fieldname='sm_fdr_gbs')]


def get_expstr_fieldnames():
    return \
"""    Available methods:
        local(
            bonferroni
            sidak
            holm
            fdr
        )
        statsmodels(
            sm_bonferroni
            sm_sidak
            holm_sidak
            sm_holm
            simes_hochberg
            hommel
            fdr_bh
            fdr_by
            fdr_tsbh
            fdr_tsbky
            fdr_gbs
        )"""

def get_exp_fieldnames():
    return cx.OrderedDict([
        (('local', 'bonferroni'), 'bonferroni'),
        (('local', 'sidak'), 'sidak'),
        (('local', 'holm'), 'holm'),
        (('local', 'fdr'), 'fdr'),
        (('statsmodels', 'bonferroni'), 'sm_bonferroni'),
        (('statsmodels', 'sidak'), 'sm_sidak'),
        (('statsmodels', 'holm-sidak'), 'holm_sidak'),
        (('statsmodels', 'holm'), 'sm_holm'),
        (('statsmodels', 'simes-hochberg'), 'simes_hochberg'),
        (('statsmodels', 'hommel'), 'hommel'),
        (('statsmodels', 'fdr_bh'), 'fdr_bh'),
        (('statsmodels', 'fdr_by'), 'fdr_by'),
        (('statsmodels', 'fdr_tsbh'), 'fdr_tsbh'),
        (('statsmodels', 'fdr_tsbky'), 'fdr_tsbky'),
        (('statsmodels', 'fdr_gbs'), 'fdr_gbs'),
    ])


def test_all():
    """Run all tests on Methods class."""
    test_init_methods()

if __name__ == '__main__':
    test_all()
