# encoding: utf-8
from __future__ import print_function, division, absolute_import

import dill


def execute(payload):
    function, args = dill.loads(payload)
    return dill.dumps(function(args))


def p_map(pool, function, data):
    """replacement for Pool.map from multiprocessing.
    The ability to pickle objects using the Python pickle / cPickle module
    limits the functions and data you can process with Pool.map.

    Here we use dill instead, which eg can pickle instance methods.

    Inspired by
    https://stackoverflow.com/questions/8804830/python-multiprocessing-pickling-error/24673524#24673524
    """

    assert isinstance(data, (list, tuple, set))
    payload = [dill.dumps((function, args)) for args in data]
    results = pool.map(execute, payload)
    return [dill.loads(result) for result in results]
