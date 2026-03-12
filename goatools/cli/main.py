"""Top-level GOATOOLS command-line dispatcher."""

from __future__ import print_function

import importlib
import sys

from goatools import __version__


SUBCOMMANDS = {
    "compare_gos": (
        "goatools.cli.compare_gos",
        "Compare two or more sets of GO IDs.",
    ),
    "fetch_associations": (
        "goatools.cli.fetch_associations",
        "Download associations from GOlr for a species.",
    ),
    "find_enrichment": (
        "goatools.cli.find_enrichment",
        "Run GO enrichment analysis.",
    ),
    "go_plot": (
        "goatools.cli.gosubdag_plot",
        "Plot GO hierarchies.",
    ),
    "map_to_slim": (
        "goatools.cli.map_to_slim",
        "Map GO terms or associations to GO slim terms.",
    ),
    "ncbi_gene_results_to_python": (
        "goatools.cli.ncbi_gene_results_to_python",
        "Convert NCBI Gene TSV results to Python modules.",
    ),
    "plot_go_term": (
        "goatools.cli.plot_go_term",
        "Plot a GO term lineage.",
    ),
    "prt_terms": (
        "goatools.cli.prt_terms",
        "Print GO terms from the DAG source.",
    ),
    "wr_hier": (
        "goatools.cli.wr_hierarchy",
        "Write a GO hierarchy report.",
    ),
    "wr_sections": (
        "goatools.cli.wr_sections",
        "Create or update GO grouping sections files.",
    ),
}


def _get_help():
    """Return the top-level CLI help text."""
    lines = [
        "Usage:",
        "  goatools <command> [<args>...]",
        "  goatools help <command>",
        "  goatools --version",
        "  goatools -h | --help",
        "",
        "Commands:",
    ]
    for name in sorted(SUBCOMMANDS):
        _, desc = SUBCOMMANDS[name]
        lines.append("  {CMD:28} {DESC}".format(CMD=name, DESC=desc))
    return "\n".join(lines)


def _run_subcommand(subcommand, args):
    """Import and run a subcommand."""
    if subcommand not in SUBCOMMANDS:
        sys.stderr.write(
            "Unknown command: {CMD}\n\n{HELP}\n".format(
                CMD=subcommand, HELP=_get_help()
            )
        )
        raise SystemExit(2)
    modname, _ = SUBCOMMANDS[subcommand]
    module = importlib.import_module(modname)
    module.main(args)


def main(argv=None):
    """Run the GOATOOLS command-line dispatcher."""
    args = sys.argv[1:] if argv is None else argv
    if args and args[0] == "--version":
        print(__version__)
        return
    if not args or args[0] in ("-h", "--help"):
        print(_get_help())
        return
    if args[0] == "help":
        if len(args) == 1:
            print(_get_help())
            return
        _run_subcommand(args[1], ["--help"])
        return
    _run_subcommand(args[0], args[1:])
