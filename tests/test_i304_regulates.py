from goatools.gosubdag.gosubdag import GoSubDag
from goatools.obo_parser import GODag


def test_i304_regulates():
    godag = GODag("go-basic.obo", optional_attrs=["relationship"])
    optional_relationships = {
        "regulates",
        "negatively_regulates",
        "positively_regulates",
    }
    anc = GoSubDag(
        ["GO:0019222"], godag, optional_relationships, prt=None
    ).rcntobj.go2ancestors["GO:0019222"]
    assert len(anc) == 4
