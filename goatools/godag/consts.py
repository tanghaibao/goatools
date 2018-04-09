"""GO-DAG constants."""

__copyright__ = "Copyright (C) 2010-2018, DV Klopfenstein, H Tang, All rights reserved."
__author__ = "DV Klopfenstein"


# pylint: disable=too-few-public-methods
class Consts(object):
    """Constants commonly used in GO-DAG operations."""

    NAMESPACE2NS = {
        'biological_process' : 'BP',
        'molecular_function' : 'MF',
        'cellular_component' : 'CC'}

    # https://owlcollab.github.io/oboformat/doc/GO.format.obo-1_4.html
    RELATIONSHIP_LIST = ['part_of', 'regulates', 'negatively_regulates', 'positively_regulates']

    top_terms = set(['GO:0008150', 'GO:0003674', 'GO:0005575']) # BP, MF, CC
    NS2GO = {'BP':'GO:0008150', 'MF':'GO:0003674', 'CC':'GO:0005575'}
    NAMESPACE2GO = {'biological_process':'GO:0008150',
                    'molecular_function':'GO:0003674',
                    'cellular_component':'GO:0005575'}

    def __init__(self):
        self.relationships = set(self.RELATIONSHIP_LIST)


# Copyright (C) 2010-2018, DV Klopfenstein, H Tang, All rights reserved.
