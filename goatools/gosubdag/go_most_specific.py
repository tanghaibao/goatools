"""Various methods to estimating if a GO term is more specific than another GO term."""

__copyright__ = "Copyright (C) 2016-2018, DV Klopfenstein, H Tang, All rights reserved."
__author__ = "DV Klopfenstein"


def get_most_specific_dcnt(goids, go2nt):
    """Get the GO ID with the lowest descendants count."""
    # go2nt_usr = {go:go2nt[go] for go in goids}
    # return min(go2nt_usr.items(), key=lambda t: t[1].dcnt)[0]
    return min(_get_go2nt(goids, go2nt), key=lambda t: t[1].dcnt)[0]

def get_most_specific_tinfo(goids, go2nt):
    """Get the GO ID with the highest GO term annotation information value."""
    # go2nt_usr = {go:go2nt[go] for go in goids}
    # return max(go2nt_usr.items(), key=lambda t: t[1].tinfo)[0]
    return max(_get_go2nt(goids, go2nt), key=lambda t: t[1].tinfo)[0]

def get_most_specific_tinfo_dcnt(goids, go2nt):
    """Get the GO ID with the highest GO term annotation information value."""
    # go2nt_usr = {go:go2nt[go] for go in goids}
    # return max(go2nt_usr.items(), key=lambda t: [t[1].tinfo, t[1].dcnt])[0]
    return max(_get_go2nt(goids, go2nt), key=lambda t: [t[1].tinfo, t[1].dcnt])[0]

def _get_go2nt(goids, go2nt_all):
    """Get user go2nt using main GO IDs, not alt IDs."""
    go_nt_list = []
    goids_seen = set()
    for goid_usr in goids:
        ntgo = go2nt_all[goid_usr]
        goid_main = ntgo.id
        if goid_main not in goids_seen:
            goids_seen.add(goid_main)
            go_nt_list.append((goid_main, ntgo))
    return go_nt_list


# Copyright (C) 2016-2018, DV Klopfenstein, H Tang, All rights reserved.
