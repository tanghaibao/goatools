"""Small lightweight utilities used frequently in GOATOOLS."""

__copyright__ = "Copyright (C) 2016-2018, DV Klopfenstein, H Tang, All rights reserved."
__author__ = "DV Klopfenstein"


def extract_kwargs(args, exp_keys, exp_elems):
    """Return user-specified keyword args in a dictionary and a set (for True/False items)."""
    arg_dict = {}    # For arguments that have values
    arg_set = set()  # For arguments that are True or False (present in set if True)
    for key, val in args.items():
        if exp_keys is not None and key in exp_keys and val:
            arg_dict[key] = val
        elif exp_elems is not None and key in exp_elems and val:
            arg_set.add(key)
    return {'dict':arg_dict, 'set':arg_set}

def get_kwargs_set(args, exp_elem2dflt):
    """Return user-specified keyword args in a dictionary and a set (for True/False items)."""
    arg_set = set()  # For arguments that are True or False (present in set if True)
    # Add user items if True
    for key, val in args.items():
        if exp_elem2dflt is not None and key in exp_elem2dflt and val:
            arg_set.add(key)
    # Add defaults if needed
    for key, dfltval in exp_elem2dflt.items():
        if dfltval and key not in arg_set:
            arg_set.add(key)
    return arg_set

def get_kwargs(args, exp_keys, exp_elems):
    """Return user-specified keyword args in a dictionary and a set (for True/False items)."""
    arg_dict = {}    # For arguments that have values
    for key, val in args.items():
        if exp_keys is not None and key in exp_keys and val:
            if isinstance(val, str):
                val = val.strip()
            arg_dict[key] = val
        elif exp_elems is not None and key in exp_elems and val:
            arg_dict[key] = True
    return arg_dict


# Copyright (C) 2016-2018, DV Klopfenstein, H Tang, All rights reserved.
