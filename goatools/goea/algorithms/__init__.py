# -*- coding: UTF-8 -*-
"""GOEA term-scoring algorithms.

    >>> from goatools.goea.algorithms import get_algorithm
    >>> get_algorithm('classic')
"""

__copyright__ = "Copyright (C) 2010-present, H Tang et al., All rights reserved."

from goatools.goea.algorithms.base import GoeaAlgorithm, GoeaAlgoResult, GoeaContext
from goatools.goea.algorithms.classic import ClassicAlgorithm

__all__ = [
    "GoeaAlgorithm",
    "GoeaAlgoResult",
    "GoeaContext",
    "ClassicAlgorithm",
    "ALGORITHMS",
    "get_algorithm",
]

ALGORITHMS = {
    ClassicAlgorithm.name: ClassicAlgorithm,
}


def get_algorithm(algorithm=None, **kws):
    """Get a GOEA algorithm from a name, an instance, or a class.

    Keyword args prefixed with the algorithm's name are passed to its
    constructor, e.g. get_algorithm('elim', elim_cutoff=.01) -> ElimAlgorithm(cutoff=.01)
    """
    if algorithm is None:
        algorithm = ClassicAlgorithm.name
    if isinstance(algorithm, GoeaAlgorithm):
        return algorithm
    if isinstance(algorithm, type) and issubclass(algorithm, GoeaAlgorithm):
        return algorithm(**_get_algo_kws(algorithm.name, kws))
    if algorithm not in ALGORITHMS:
        raise ValueError(
            "UNKNOWN GOEA ALGORITHM({A}); EXPECTED ONE OF: {E}".format(
                A=algorithm, E=" ".join(sorted(ALGORITHMS))
            )
        )
    cls = ALGORITHMS[algorithm]
    return cls(**_get_algo_kws(cls.name, kws))


def _get_algo_kws(name, kws):
    """Pull '<name>_<arg>' keywords out of a GOEnrichmentStudy's kwargs."""
    pre = "{N}_".format(N=name)
    return {k[len(pre):]: v for k, v in kws.items() if k.startswith(pre)}
