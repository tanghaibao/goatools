#!/usr/bin/env python
"""Test notebook code"""

from goatools.base import get_godag
from goatools.gosubdag.gosubdag import GoSubDag
from goatools.gosubdag.plot.gosubdag_plot import GoSubDagPlot

def test_nb():
    """Test notebook code"""
    godag = get_godag("go-basic.obo", optional_attrs={'relationship'})
    go_leafs = set(o.item_id for o in godag.values() if not o.children)
    virion = 'GO:0019012'
    gosubdag_r0 = GoSubDag(go_leafs, godag)
    nt_virion = gosubdag_r0.go2nt[virion]
    print(nt_virion)
    print('THE VALUE OF dcnt IS: {dcnt}'.format(dcnt=nt_virion.dcnt))
    gosubdag_r1 = GoSubDag(go_leafs, godag, relationships=True)
    nt_virion = gosubdag_r1.go2nt[virion]
    print(nt_virion)
    print('THE VALUE OF dcnt IS: {dcnt}'.format(dcnt=nt_virion.dcnt))
    gosubdag_partof = GoSubDag(go_leafs, godag, relationships={'part_of'})
    nt_virion = gosubdag_partof.go2nt[virion]
    print(nt_virion)
    print('THE VALUE OF dcnt IS: {dcnt}'.format(dcnt=nt_virion.dcnt))
    virion_descendants = gosubdag_partof.rcntobj.go2descendants[virion]
    print('{N} descendants of virion were found'.format(N=len(virion_descendants)))

    # Limit plot of descendants to get a smaller plot
    virion_capsid_fiber = {'GO:0098033', 'GO:0098032'}
    gosubdag_partof.prt_goids(virion_capsid_fiber,
                              '{NS} {GO} dcnt({dcnt}) D-{depth:02} {GO_name}')

    # Limit plot size by choosing just two virion descendants
    # Get a subset containing only a couple virion descendants and their ancestors
    pltdag = GoSubDag(virion_capsid_fiber, godag, relationships={'part_of'})
    pltobj = GoSubDagPlot(pltdag)
    pltobj.plt_dag('virion_capsid_fiber.png')

if __name__ == '__main__':
    test_nb()
