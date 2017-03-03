"""Runs stats describe on data. Prints data in a markdown table."""

__copyright__ = "Copyright (C) 2016-2017, DV Klopfenstein, H Tang. All rights reserved."
__author__ = "DV Klopfenstein"

import sys
import math
import numpy as np
import scipy.stats as stats

class StatsDescribe(object):
    """Describe summary statistics for a list of numbers."""

    fmt = "{name:10} | {qty:13} | {range:20} | {25th percentile:>15} | " \
          "{median:>8} | {75th percentile:>15} | {mean:>8} | {stddev:>}\n"

    def __init__(self, desc, fmtstr="{:>8.2e}"):
        self.desc = desc
        self.fmtstr = fmtstr

    def prt_hdr(self, prt=sys.stdout):
        """Print stats header in markdown style."""
        hdr = "name       | # {items:11} | range of {items:11} | 25th percentile | " \
              "  median | 75th percentile |     mean | stddev\n".format(items=self.desc)
        div = "-----------|---------------|----------------------|" \
              "-----------------|----------|-----------------|----------|-------\n"
        prt.write("\n{DESC}\n".format(DESC=self.desc))
        prt.write(hdr)
        prt.write(div)

    def prt_data(self, name, vals, prt=sys.stdout):
        """Print stats data in markdown style."""
        fld2val = self._init_fld2val(name, vals)
        prt.write(self.fmt.format(**fld2val))
        return fld2val

    def _init_fld2val(self, name, vals):
        """Describe summary statistics for a list of numbers."""
        if vals:
            return self._init_fld2val_stats(name, vals)
        return self._init_fld2val_null(name)

    def _init_fld2val_stats(self, name, vals):
        """Return statistics on values."""
        vals_stats = stats.describe(vals)
        stddev = math.sqrt(vals_stats[3]) # stats variance
        p25 = np.percentile(vals, 25)
        p50 = np.percentile(vals, 50) # median
        p75 = np.percentile(vals, 75)
        fld2val = {
            'name':name,
            'qty'.format(ITEMS=self.desc):vals_stats[0], # stats nobs
            'range':self._get_str_range(vals_stats),
            '25th percentile':p25,
            'median':p50,
            '75th percentile':p75,
            'mean':vals_stats[2], # stats mean
            'stddev':stddev}
        fmtflds = set(['25th percentile', 'median', '75th percentile', 'mean', 'stddev'])
        mkint = "," in self.fmtstr
        for key, val in fld2val.items():
            if key in fmtflds:
                if mkint:
                    val = int(round(val))
                fld2val[key] = self.fmtstr.format(val)
        return fld2val

    def _init_fld2val_null(self, name):
        """Return empty fields if there are no values."""
        return {
            'name':name,
            'qty'.format(ITEMS=self.desc):0,
            'range':"",
            '25th percentile':"",
            'median':"",
            '75th percentile':"",
            'mean':"",
            'stddev':""}

    def _get_str_range(self, vals_stats):
        """Return a string containing the range of values."""
        minmax = vals_stats[1] # stats minmax
        minval = self.fmtstr.format(minmax[0])
        maxval = self.fmtstr.format(minmax[1])
        return '{A} to {B:6>}'.format(A=minval, B=maxval)

# Copyright (C) 2016-2017, DV Klopfenstein, H Tang. All rights reserved.
