"""GO-DAG constants."""

__copyright__ = "Copyright (C) 2010-2019, DV Klopfenstein, H Tang, All rights reserved."
__author__ = "DV Klopfenstein"


NAMESPACE2NS = {
    'biological_process' : 'BP',
    'molecular_function' : 'MF',
    'cellular_component' : 'CC'}

NS2NAMESPACE = {
    'BP': 'biological_process',
    'MF': 'molecular_function',
    'CC': 'cellular_component'}

# https://owlcollab.github.io/oboformat/doc/GO.format.obo-1_4.html
RELATIONSHIP_LIST = ['part_of', 'regulates', 'negatively_regulates', 'positively_regulates']
RELATIONSHIP_SET = set(RELATIONSHIP_LIST)

TOP_TERMS = set(['GO:0008150', 'GO:0003674', 'GO:0005575']) # BP, MF, CC
NS2GO = {'BP':'GO:0008150', 'MF':'GO:0003674', 'CC':'GO:0005575'}
NAMESPACE2GO = {'biological_process':'GO:0008150',
                'molecular_function':'GO:0003674',
                'cellular_component':'GO:0005575'}

def chk_relationships(relationships):
    """Check if user-provided relationships are seen in the GO DAG"""
    assert relationships.issubset(RELATIONSHIP_SET), 'RELATIONSHIP({R}) NOT IN: {Rs}'.format(
        R=relationships.difference(RELATIONSHIP_SET), Rs=RELATIONSHIP_SET)


# Copyright (C) 2010-2019, DV Klopfenstein, H Tang, All rights reserved.
