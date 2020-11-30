"""Set a random seed."""

__copyright__ = "Copyright (C) 2015-present, DV Klopfenstein, Haibao Tang. All rights reserved."
__author__ = "DV Klopfenstein"

import sys
import numpy as np


class RandomSeed32:
    """Manage a random seed."""

    def __init__(self, seed=None):
        self.max_val = int("0xffffffff", 0)
        self.seed = self._init_seed(seed)

    def _init_seed(self, seed):
        """Use user-specified random seed or randomly pick one."""
        # None -> Create seed
        clsname = self.__class__.__name__
        if seed is None:
            seed = np.random.randint(self.max_val)
            sys.stdout.write("  RANDOM SEED({S}) CREATED BY {CLASS}\n".format(
                S=self.get_seed_str(seed), CLASS=clsname))
        # Set seed
        np.random.seed(seed)
        sys.stdout.write("  RANDOM SEED({S}) STORED IN {CLASS}\n".format(
            S=self.get_seed_str(seed), CLASS=clsname))
        return seed

    def prt(self, prt=sys.stdout):
        """Print object's random seed."""
        self.prt_seed(self.seed, prt)

    def prt_seed(self, seed, prt=sys.stdout):
        """Print given random seed."""
        prt.write("  RANDOM SEED = {SEED}\n".format(SEED=self.get_seed_str(seed)))

    @staticmethod
    def get_seed_str(seed):
        """Print given random seed."""
        return "0x{S:08x} = {S:,}".format(S=seed)

    def get_seed_hexstr(self):
        """Return a string representing the random seed in hex format"""
        return '0x{S:08x}'.format(S=self.seed)


# Copyright (C) 2015-present, DV Klopfenstein, Haibao Tang. All rights reserved.
