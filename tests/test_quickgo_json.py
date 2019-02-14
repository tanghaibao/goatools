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

    print('\nREAD JSON RESULTS')
    jsonstr = obj.read()
    # print(jsonstr)
    print(type(jsonstr))

    jsonout = json.loads(jsonstr)
    print(jsonout.keys())
    print('numberOfHits: {V}'.format(V=jsonout['numberOfHits']))
    print('pageInfo: {V}'.format(V=jsonout['pageInfo']))

    # PRINT THE GO TERM
    print('\nresults:')
    assert len(jsonout['results']) == 1
    go_dict = jsonout['results'][0]
    print('id: {GO}'.format(GO=go_dict['id']))
    print('{N} children'.format(N=len(go_dict['children'])))
    print('DEFINITION: {DEFN}'.format(DEFN=go_dict['definition']))
    print(jsonout['results'][0].keys())


if __name__ == '__main__':
    test_quickgo()

# Copyright (C) 2015-2019, DV Klopfenstein, H Tang, All rights reserved.
