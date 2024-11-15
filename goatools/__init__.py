from datetime import datetime
from importlib.metadata import version, PackageNotFoundError

__author__ = ("Haibao Tang", "DV Klopfenstein")
__copyright__ = (
    f"Copyright (C) 2009-{datetime.now().year}, Haibao Tang, DV Klopfenstein"
)
__email__ = "tanghaibao@gmail.com"
__license__ = "BSD"
__status__ = "Development"

try:
    VERSION = version(__name__)
except PackageNotFoundError:  # pragma: no cover
    try:
        from .version import version as VERSION  # noqa
    except ImportError as exc:  # pragma: no cover
        raise ImportError(
            "Failed to find (autogenerated) version.py. "
            "This might be because you are installing from GitHub's tarballs, "
            "use the PyPI ones."
        ) from exc
__version__ = VERSION
