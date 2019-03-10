#!/usr/bin/env python
"""Download a oboxml file from QuickGO and read it."""

from __future__ import print_function

__copyright__ = "Copyright (C) 2015-2019, DV Klopfenstein, H Tang, All rights reserved."

import pytest

@pytest.mark.skip
def test_quickgo():
    """Download a oboxml file from QuickGO and read it."""
    from urllib import urlopen
    import xmltodict

    # NOTE: Use json instead of xml. See test_quickgo_json.py
    oboxml = urlopen('http://www.ebi.ac.uk/QuickGO/GTerm?id=GO:0003723&format=oboxml')
    print(oboxml)

    # __doc__ __init__ __iter__ __module__ __repr__
    # close code fileno fp getcode geturl headers info next read readline readlines url
    print(' '.join(sorted(dir(oboxml))))

    # Code
    print(oboxml.code)
    assert oboxml.code == 200

    # Headers
    # Example:
    #     Server: Apache/2.2.15 (Red Hat)
    #     Content-Type: text/html; charset=UTF-8
    #     Strict-Transport-Security: max-age=0
    #     Date: Thu, 14 Feb 2019 01:59:46 GMT
    #     Accept-Ranges: bytes
    #     ETag: "fc747f73f-181c-57bb822b46580"
    #     Connection: close
    #     Last-Modified: Wed, 28 Nov 2018 11:47:50 GMT
    #     Content-Length: 6172
    print(oboxml.headers)

    # <bound method addinfourl.info of <addinfourl at 7696577623752 whose
    # fp = <socket._fileobject object at 0x6ffffc601d0>>>
    print(oboxml.info)

    print('\nREAD')
    htmlstr = oboxml.read()
    print(htmlstr)
    print(type(htmlstr))
    return
    obo_dict_parse = xmltodict.parse(htmlstr)
    print(obo_dict_parse)


if __name__ == '__main__':
    test_quickgo()

# Copyright (C) 2015-2019, DV Klopfenstein, H Tang, All rights reserved.
