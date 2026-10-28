"""Sphinx configuration for the GOATOOLS documentation site.

Build locally with `make -C sphinx html`; the site is in sphinx/_build/html.
"""

import re
import shutil
import sys
from pathlib import Path
from typing import List

DOCS = Path(__file__).resolve().parent
ROOT = DOCS.parent
sys.path.insert(0, str(ROOT))


def _github_alerts_to_admonitions(text: str) -> str:
    """Turn GitHub's "> [!TIP]" callouts into MyST admonitions, which Sphinx renders."""
    lines, out, i = text.split("\n"), [], 0
    while i < len(lines):
        mtch = re.match(r"> \[!(NOTE|TIP|IMPORTANT|WARNING|CAUTION)\]\s*$", lines[i])
        if not mtch:
            out.append(lines[i])
            i += 1
            continue
        out.append("```{" + mtch.group(1).lower() + "}")
        i += 1
        while i < len(lines) and lines[i].startswith(">"):
            out.append(lines[i][1:].removeprefix(" "))
            i += 1
        out.append("```")
    return "\n".join(out)


# The site reuses content kept elsewhere in the repository. Sphinx only reads
# documents inside this directory, so copy it in on every build (gitignored).
# Copies, not symlinks: Sphinx resolves symlinked link targets back to the
# repository root, outside this directory, and README.md's links would break.
for name in ("notebooks", "docs", "README.md"):
    src, dst = ROOT / name, DOCS / name
    if dst.is_symlink():
        dst.unlink()
    elif dst.is_dir():
        shutil.rmtree(dst)
    if src.is_dir():
        shutil.copytree(src, dst, ignore=shutil.ignore_patterns(".ipynb_checkpoints"))
    else:
        dst.write_text(_github_alerts_to_admonitions(src.read_text()))

from goatools import __version__  # noqa: E402

project = "GOATOOLS"
author = "Haibao Tang et al."
copyright = "2010-present, Haibao Tang et al."
version = release = __version__

extensions = [
    "myst_nb",
    "sphinx.ext.autodoc",
    "sphinx.ext.napoleon",
    "sphinx.ext.viewcode",
    "sphinx.ext.intersphinx",
]

# Notebooks are rendered with their saved outputs; nothing is executed
nb_execution_mode = "off"
myst_enable_extensions = ["colon_fence"]
myst_heading_anchors = 4

# README.md is part of index.md, not a page of its own
exclude_patterns = ["_build", "README.md"]

autodoc_default_options = {"members": True, "show-inheritance": True}
autodoc_member_order = "bysource"
# Show defaults as written (prt=sys.stdout), not as their repr
autodoc_preserve_defaults = True
intersphinx_mapping = {"python": ("https://docs.python.org/3", None)}

# Literal blocks without a language (usage text, ASCII tables) stay plain
highlight_language = "none"

html_theme = "furo"
html_title = "GOATOOLS"


def _cli_usage_verbatim(app, what: str, name: str, obj, options, lines: List[str]) -> None:
    """Show goatools.cli module docstrings verbatim.

    They are command-line usage text, parsed as such by DocOptParse and argparse,
    not reStructuredText.
    """
    if what == "module" and name.startswith("goatools.cli.") and lines:
        lines[:] = ["::", ""] + ["    " + line for line in lines]


def setup(app) -> None:
    app.connect("autodoc-process-docstring", _cli_usage_verbatim)
