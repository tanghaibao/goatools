"""Functions for shortening GO names."""

__copyright__ = "Copyright (C) 2016-2017, DV Klopfenstein, H Tang, All rights reserved."
__author__ = "DV Klopfenstein"


greek2uni = {
    "alpha": "α",
    "beta": "β",
    "gamma": "γ",
    "delta": "δ",
}
greek2tex = {
    "alpha": r"$\alpha$",
    "beta": r"$\beta$",
    "gamma": r"$\gamma$",
    "delta": r"$\delta$",
}
keep = [
    "defense response to protozoan",
    "defense response to bacterium",
    "cellular response to interferon-beta",
    "defense response to virus",
    "response to interferon-gamma",
    "innate immune response",
    "inflammatory response",
    "response to virus",
    "immune response",
]


def get_short_plot_name(goobj):
    """Shorten some GO names so plots are smaller."""
    name = goobj.name
    if _keep_this(name):
        return replace_greek(name)
    name = name.replace(
        "cellular response to chemical stimulus", "cellular rsp. to chemical stim."
    )
    depth = goobj.depth
    if depth > 1:
        name = name.replace("regulation of", "reg. of")
        name = name.replace("positive reg", "+reg")
        name = name.replace("negative reg", "-reg")
        name = name.replace("involved in", "in")
    if depth > 2:
        name = name.replace("antigen processing and presentation", "a.p.p")
        name = name.replace("MHC class I", "MHC-I")
        if depth == 4:
            if goobj.id == "GO:0002460":
                before = " ".join(
                    [
                        "adaptive immune response based on somatic recombination of",
                        "immune receptors built from immunoglobulin superfamily domains",
                    ]
                )
                name = name.replace(
                    before,
                    "rsp. based on somatic recombination of Ig immune receptors",
                )
        if depth > 3:
            name = name.replace("signaling pathway", "sig. pw.")
            name = name.replace("response", "rsp.")
            name = name.replace("immunoglobulin superfamily domains", "Ig domains")
            name = name.replace("immunoglobulin", "Ig")
        if depth > 4:
            name = name.replace("production", "prod.")
        if depth == 6 or depth == 5:
            name = name.replace("tumor necrosis factor", "TNF")
    name = replace_greek(name)
    return name


def shorten_go_name_ptbl1(name):
    """Shorten GO name for tables in paper."""
    if _keep_this(name):
        return name
    name = name.replace("negative", "neg.")
    name = name.replace("positive", "pos.")
    name = name.replace("response", "rsp.")
    name = name.replace("regulation", "reg.")
    name = name.replace("antigen processing and presentation", "app.")
    return name


def shorten_go_name_ptbl3(name, dcnt):
    """Shorten GO description for Table 3 in manuscript."""
    if _keep_this(name):
        return name
    name = name.replace(
        "positive regulation of immune system process",
        "+ reg. of immune sys. process",
    )
    name = name.replace(
        "positive regulation of immune response", "+ reg. of immune response"
    )
    name = name.replace(
        "positive regulation of cytokine production",
        "+ reg. of cytokine production",
    )
    if dcnt < 40:
        name = name.replace("antigen processing and presentation", "a.p.p.")
    if dcnt < 10:
        name = name.replace("negative", "-")
        name = name.replace("positive", "+")
        name = name.replace("tumor necrosis factor production", "TNF production")
    if dcnt < 4:
        name = name.replace("regulation", "reg.")
        name = name.replace("exogenous ", "")
        name = name.replace(" via ", " w/")
        name = name.replace("T cell mediated cytotoxicity", "cytotoxicity via T cell")
    name = name.replace("involved in", "in")
    name = name.replace("-positive", "+")
    return name


def replace_greek(name):
    """Replace text representing greek letters with greek letters."""
    name = name.replace("gamma-delta", "gammadelta")
    name = name.replace("interleukin-1 beta", "interleukin-1beta")
    for greek_txt, uni in greek2uni.items():
        if greek_txt in name:
            name = name.replace(greek_txt, uni)
    return name


def replace_greek_tex(name):
    """Replace text representing greek letters with greek letters."""
    name = name.replace("gamma-delta", "gammadelta")
    name = name.replace("interleukin-1 beta", "interleukin-1beta")
    for greek_txt, tex in greek2tex.items():
        if greek_txt in name:
            name = name.replace(greek_txt, tex)
    return name


def shorten_go_name_all(name):
    """Shorten GO name for tables in paper, supplemental materials, and plots."""
    name = replace_greek(name)
    name = name.replace("MHC class I", "MHC-I")
    return name


def _keep_this(name):
    """Return True if there are to be no modifications to name."""
    for keep_name in keep:
        if name == keep_name:
            return True
    return False


# Copyright (C) 2016-2017, DV Klopfenstein, H Tang, All rights reserved.
