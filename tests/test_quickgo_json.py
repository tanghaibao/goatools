#!/usr/bin/env python
"""Download a obj file from QuickGO and read it."""

from __future__ import print_function

__copyright__ = "Copyright (C) 2015-2019, DV Klopfenstein, H Tang, All rights reserved."


from urllib import urlopen
import json

def test_quickgo():
    """Download a obj file from QuickGO and read it."""
    obj = urlopen('https://www.ebi.ac.uk/QuickGO/services/ontology/go/terms/GO:0003723/complete')
    print(obj)

    # __doc__ __init__ __iter__ __module__ __repr__
    # close code fileno fp getcode geturl headers info next read readline readlines url
    print(' '.join(sorted(dir(obj))))

    # Code
    print(obj.code)
    assert obj.code == 200

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
    print(obj.headers)

    # <bound method addinfourl.info of <addinfourl at 7696577623752 whose
    # fp = <socket._fileobject object at 0x6ffffc601d0>>>
    print(obj.info)

    print('\nREAD')
    jsonstr = obj.read()
    # print(jsonstr)
    print(type(jsonstr))

    go_dict = json.loads(jsonstr)
    print(go_dict.keys())
    print('numberOfHits: {V}'.format(V=go_dict['numberOfHits']))
    print('pageInfo: {V}'.format(V=go_dict['pageInfo']))
    print('results:')
    assert len(go_dict['results']) == 1
    for key, val in sorted(go_dict['results'][0].items()):
        print('    {N:3} {K}'.format(K=key, N=len(val)))


if __name__ == '__main__':
    test_quickgo()

# Copyright (C) 2015-2019, DV Klopfenstein, H Tang, All rights reserved.
