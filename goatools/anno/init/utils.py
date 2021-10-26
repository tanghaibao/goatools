"""Utilities for combining namedtuples."""

__copyright__ = "Copyright (C) 2016-2019, DV Klopfenstein, H Tang. All rights reserved."
__author__ = "DV Klopfenstein"

from datetime import date


def get_date_yyyymmdd(yyyymmdd):
    """Return datetime.date given string."""
    return date(int(yyyymmdd[:4]), int(yyyymmdd[4:6], base=10), int(yyyymmdd[6:], base=10))


# Copyright (C) 2016-2019, DV Klopfenstein, H Tang. All rights reserved.
