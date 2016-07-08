"""Tasks for GOTerms in obo dag."""

def get_all_parents(go_objs):
    """Return a set containing all GO Term parents of multiple GOTerm objects."""
    go_parents = set()
    for go_obj in go_objs:
        go_parents |= go_obj.get_all_parents()
    return go_parents

