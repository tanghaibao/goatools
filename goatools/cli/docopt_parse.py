"""Run docopt in GOATOOLS."""

from __future__ import print_function

__copyright__ = "Copyright (C) 2016-2019, DV Klopfenstein, H Tang, All rights reserved."
__author__ = "DV Klopfenstein"

import sys
import re
from docopt import docopt
from goatools.gosubdag.utils import get_kwargs


#pylint: disable=too-few-public-methods
class DocOptParse(object):
    """Run docopt in GOATOOLS."""

    def __init__(self, doc, exp_keys, exp_elems):
        self.doc = doc  # doc string
        self.exp_keys = exp_keys if exp_keys else set()     # Expected dictionary keys
        self.exp_elems = exp_elems if exp_elems else set()  # Expected set elements (True/False)

    def get_docargs(self, args=None, prt=None, **kws):
        """Pare down docopt. Return a minimal dictionary and a set containing runtime arg values."""
        arg_kws = self._get_docargs(args, prt)
        if 'intvals' in kws:
            self._set_intvals(arg_kws, kws['intvals'])
        return arg_kws

    def _get_docargs(self, args_user, prt):
        """Pare down docopt. Return a minimal dictionary and a set containing runtime arg values."""
        if prt is not None:
            print("DocOptParse BEFORE docopt: {}".format(args_user))
        docargs = docopt(self.doc, args_user)
        if prt is not None:
            print("DocOptParse AFTER  docopt: {}".format(docargs))
        kwargs_doc = {re.sub(r'^-{1,2}', '', k):v for k, v in docargs.items()}
        self._chk_docopt_kws(kwargs_doc, args_user)
        kwargs_usr = get_kwargs(kwargs_doc, self.exp_keys, self.exp_elems)
        if prt is not None:
            print("DocOptParse AFTER  pared: {}".format(kwargs_usr))
        for key in ['taxid']:
            if key in kwargs_usr:
                kwargs_usr[key] = int(kwargs_usr[key])
        if prt is not None:
            print("DocOptParse AFTER  edited/checked: {}".format(kwargs_usr))
        return kwargs_usr

    @staticmethod
    def _set_intvals(kws, keys):
        """Convert keyword values to int."""
        for key in keys:
            if key in kws:
                kws[key] = int(kws[key])

    def _chk_docopt_exit(self, args, exp_letters):
        """Check if docopt exit was for an unknown argument."""
        if args is None:
            args = sys.argv[1:]
        keys_all = self.exp_keys.union(self.exp_elems)
        if exp_letters:
            keys_all |= exp_letters
        unknown_args = self._chk_docunknown(args, keys_all)
        if unknown_args:
            raise RuntimeError("{USAGE}\n    **FATAL: UNKNOWN ARGS: {UNK}".format(
                USAGE=self.doc, UNK=" ".join(unknown_args)))

    def _chk_docopt_kws(self, docdict, exp):
        """Check for common user errors when running from the command-line."""
        for key, val in docdict.items():
            if isinstance(val, str):
                assert '=' not in val, self._err("'=' FOUND IN VALUE", key, val, exp)
            elif key != 'help' and key not in self.exp_keys and key not in self.exp_elems:
                raise RuntimeError(self._err("UNKNOWN KEY", key, val, exp))

    def _err(self, msg, key, val, exp):
        return "{DOC}\n{MSG}: KEY({K}) VAL({V}): {EXP}".format(
            DOC=self.doc, MSG=msg, K=key, V=val, EXP=" ".join(exp))

    @staticmethod
    def _chk_docunknown(args, exp):
        """Return any unknown args."""
        unknown = []
        for arg in args:
            if arg[:2] == '--':
                val = arg[2:]
                if val not in exp:
                    unknown.append(arg)
            elif arg[:1] == '-':
                val = arg[1:]
                if val not in exp:
                    unknown.append(arg)
        if '-h' in unknown or '--help' in unknown:
            return []
        return unknown


# Copyright (C) 2016-2019, DV Klopfenstein, H Tang, All rights reserved.
