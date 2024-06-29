import pytest

from goatools.gosubdag.plot.go_name_shorten import get_short_plot_name


@pytest.mark.parametrize(
    "inp, outp",
    [
        ("inflammatory response", "inflammatory response"),
        ("response to virus", "response to virus"),
        ("immune response", "immune response"),
        ("cellular response to chemical stimulus", "cellular rsp. to chemical stim."),
        ("regulation of", "reg. of"),
        ("positive reg", "+reg"),
        ("negative reg", "-reg"),
        ("involved in", "in"),
        ("antigen processing and presentation", "a.p.p"),
        ("MHC class I", "MHC-I"),
        ("signaling pathway", "sig. pw."),
        ("response", "rsp."),
        ("immunoglobulin superfamily domains", "Ig domains"),
        ("immunoglobulin", "Ig"),
        ("production", "prod."),
        ("tumor necrosis factor", "TNF"),
        ("alpha beta T cell activation", "α β T cell activation"),
    ],
)
def test_get_short_plot_name(inp: str, outp: str):
    goobj = type("goobj", (object,), {"name": inp, "depth": 6, "id": "GO:0000000"})
    assert get_short_plot_name(goobj) == outp
