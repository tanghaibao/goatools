"""Test issue #339: _chk_goids_notfound should honor prt=None.

https://github.com/tanghaibao/goatools/issues/339

Previously _chk_goids_notfound wrote to sys.stderr unconditionally, so
GOEnrichmentStudy(..., log=None) would still emit the "GO IDs NOT FOUND IN
ASSOCIATION" message. The fix routes the message through the prt argument.
"""

import io
from contextlib import redirect_stderr, redirect_stdout

from goatools.anno.update_association import _chk_goids_notfound


def test_i339_chk_goids_notfound_silent_when_prt_none():
    """_chk_goids_notfound(prt=None) must not write to stdout or stderr."""
    goids_assoc = {"GO:0000001", "GO:9999999"}  # one missing
    goids_avail = {"GO:0000001"}

    out, err = io.StringIO(), io.StringIO()
    with redirect_stdout(out), redirect_stderr(err):
        _chk_goids_notfound(goids_assoc, goids_avail, prt=None)

    assert out.getvalue() == "", f"unexpected stdout: {out.getvalue()!r}"
    assert err.getvalue() == "", f"unexpected stderr: {err.getvalue()!r}"


def test_i339_chk_goids_notfound_writes_to_prt():
    """When a prt stream is supplied, the missing-GO message is written there."""
    goids_assoc = {"GO:0000001", "GO:9999999"}
    goids_avail = {"GO:0000001"}

    prt = io.StringIO()
    _chk_goids_notfound(goids_assoc, goids_avail, prt=prt)

    assert "GO:9999999" in prt.getvalue()


if __name__ == "__main__":
    test_i339_chk_goids_notfound_silent_when_prt_none()
    test_i339_chk_goids_notfound_writes_to_prt()
