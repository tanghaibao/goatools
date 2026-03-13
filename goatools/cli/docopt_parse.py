"""Internal CLI parsing helpers used by GOATOOLS."""

from __future__ import print_function

__copyright__ = "Copyright (C) 2016-present, DV Klopfenstein, H Tang, All rights reserved."
__author__ = "DV Klopfenstein"

import argparse
import shlex
import sys

from goatools.godag.consts import NS2GO
from goatools.gosubdag.utils import get_kwargs


class _ArgumentParser(argparse.ArgumentParser):
    """Argparse wrapper which raises runtime errors with the legacy help text."""

    def __init__(self, doc, prog):
        super().__init__(prog=prog, add_help=False, allow_abbrev=False)
        self._doc = doc

    def error(self, message):
        raise RuntimeError("{DOC}\n{MSG}\n".format(DOC=self._doc, MSG=message))


# pylint: disable=too-few-public-methods
class DocOptParse(object):
    """Parse the handful of legacy CLI modules which previously used docopt."""

    def __init__(self, doc, exp_keys, exp_elems):
        self.doc = doc
        self.exp_keys = exp_keys if exp_keys else set()
        self.exp_elems = exp_elems if exp_elems else set()
        self.cmd = self._get_command() if doc is not None else None

    def get_docargs(self, args=None, prt=None, **kws):
        """Return a pared-down runtime kwargs dictionary."""
        args_user = self._normalize_args(args)
        if prt is not None:
            print("DocOptParse BEFORE parse: {}".format(args_user))
        kwargs_usr = self._get_docargs(args_user)
        if "intvals" in kws:
            self._set_intvals(kwargs_usr, kws["intvals"])
        if prt is not None:
            print("DocOptParse AFTER  pared: {}".format(kwargs_usr))
        return kwargs_usr

    def _get_docargs(self, args_user):
        """Parse command-line arguments for a supported CLI."""
        if "-h" in args_user or "--help" in args_user:
            print(self.doc)
            raise SystemExit(0)
        argdict = self._parse_args(args_user)
        kwargs_usr = get_kwargs(argdict, self.exp_keys, self.exp_elems)
        if "taxid" in kwargs_usr:
            kwargs_usr["taxid"] = int(kwargs_usr["taxid"])
        return kwargs_usr

    def _parse_args(self, args_user):
        """Dispatch parsing to the parser for the current command."""
        cmd2parser = {
            "compare_gos": self._parse_compare_gos,
            "go_plot": self._parse_go_plot,
            "prt_terms": self._parse_prt_terms,
            "wr_hier": self._parse_wr_hier,
            "wr_sections": self._parse_wr_sections,
        }
        if self.cmd not in cmd2parser:
            raise RuntimeError("Unsupported CLI parser: {}".format(self.cmd))
        return cmd2parser[self.cmd](args_user)

    @staticmethod
    def _set_intvals(kws, keys):
        """Convert keyword values to int."""
        for key in keys:
            if key in kws:
                kws[key] = int(kws[key])

    def _normalize_args(self, args_user):
        """Split compact option/value strings and drop program/subcommand prefixes."""
        args = self._split_args(sys.argv[1:] if args_user is None else args_user)
        while args and args[0] == "goatools":
            args = args[1:]
        if self.cmd and args and args[0] == self.cmd:
            args = args[1:]
        return args

    @staticmethod
    def _split_args(args_user):
        """Split compact option/value strings while leaving positional args untouched."""
        if isinstance(args_user, str):
            return shlex.split(args_user)
        args = []
        for arg in args_user:
            if isinstance(arg, str) and arg.startswith("-") and any(ch.isspace() for ch in arg):
                args.extend(shlex.split(arg))
            else:
                args.append(arg)
        return args

    def _get_command(self):
        """Infer the subcommand name from the usage section."""
        for line in self.doc.splitlines():
            line = line.strip()
            if line.startswith("goatools "):
                parts = line.split()
                if len(parts) > 1:
                    return parts[1]
        raise RuntimeError("Unable to determine CLI command from usage text")

    def _parse_compare_gos(self, args):
        parser = _ArgumentParser(self.doc, "goatools compare_gos")
        parser.add_argument("GO_FILE", nargs="+")
        parser.add_argument("-s", "--sections")
        parser.add_argument("-S")
        parser.add_argument("-o", "--ofile")
        parser.add_argument("--xlsx")
        parser.add_argument("-v", "--verbose", action="store_true")
        parser.add_argument("--obo", default="go-basic.obo")
        parser.add_argument("--slims", default="goslim_generic.obo")
        parser.add_argument("--gaf")
        parser.add_argument("--gene2go")
        parser.add_argument("--taxid", type=int)
        parsed = vars(parser.parse_args(args))
        # nargs="+" ensures >=1; this enforces the domain rule of >=2 files
        if len(parsed["GO_FILE"]) < 2:
            parser.error("at least two GO_FILE arguments are required")
        return parsed

    def _parse_go_plot(self, args):
        parser = _ArgumentParser(self.doc, "goatools go_plot")
        parser.add_argument("GO", nargs="*")
        parser.add_argument("-i", "--go_file")
        parser.add_argument("-o", "--outfile", default="go_plot.png")
        parser.add_argument("-r", "--relationship", action="store_true")
        parser.add_argument("--relationships")
        parser.add_argument("-s", "--sections")
        parser.add_argument("-S")
        parser.add_argument("--gpad")
        parser.add_argument("--gaf")
        parser.add_argument("--id2gos")
        parser.add_argument("--gene2go")
        parser.add_argument("--taxid", type=int)
        parser.add_argument("--obo", default="go-basic.obo")
        parser.add_argument("-t", "--title")
        parser.add_argument("-p", "--parentcnt", action="store_true")
        parser.add_argument("-c", "--childcnt", action="store_true")
        parser.add_argument("--shorten", action="store_true")
        parser.add_argument("--no_ldr", action="store_true")
        parser.add_argument("--mark_alt_id", action="store_true")
        parser.add_argument("--draw-children", dest="draw-children", action="store_true")
        parser.add_argument("--go_aliases")
        parser.add_argument("--go_color_file")
        parser.add_argument("--rankdir")
        parser.add_argument("--norel", action="store_true")
        return vars(parser.parse_args(args))

    def _parse_prt_terms(self, args):
        parser = _ArgumentParser(self.doc, "goatools prt_terms")
        parser.add_argument("items", nargs="*")
        parser.add_argument("-i", "--ifile", default="sections_in.txt")
        parser.add_argument("-n", "--name")
        parser.add_argument("--obo", default="go-basic.obo")
        parsed = vars(parser.parse_args(args))
        gos = []
        go_file = None
        for item in parsed.pop("items"):
            if item.startswith("GO:") or item in NS2GO:
                gos.append(item)
            elif go_file is None:
                go_file = item
            else:
                raise RuntimeError("{DOC}\nunrecognized arguments: {ARG}\n".format(
                    DOC=self.doc, ARG=item))
        if gos:
            parsed["GO"] = gos
        if go_file is not None:
            parsed["GO_FILE"] = go_file
        return parsed

    def _parse_wr_hier(self, args):
        parser = _ArgumentParser(self.doc, "goatools wr_hier")
        parser.add_argument("GO", nargs="*")
        parser.add_argument("-i")
        parser.add_argument("-o")
        parser.add_argument("-f", action="store_true")
        parser.add_argument("--up", action="store_true")
        parser.add_argument("--dag", default="go-basic.obo")
        parser.add_argument("--gaf")
        parser.add_argument("--gene2go")
        parser.add_argument("--taxid", type=int)
        parser.add_argument("--no_indent", action="store_true")
        parser.add_argument("--max_indent", type=int)
        parser.add_argument("--concise", action="store_true")
        parser.add_argument("--dash_len", type=int, default=6)
        parser.add_argument("--item_marks")
        parser.add_argument("--include_only")
        parser.add_argument("-r", "--relationship", action="store_true")
        return vars(parser.parse_args(args))

    def _parse_wr_sections(self, args):
        parser = _ArgumentParser(self.doc, "goatools wr_sections")
        parser.add_argument("GO_FILE", nargs="?")
        parser.add_argument("-i", "--ifile", default="sections_in.txt")
        parser.add_argument("-o", "--ofile", default="sections.txt")
        parser.add_argument("--txt", default="grouped_gos.txt")
        parser.add_argument("--py")
        parser.add_argument("--xlsx")
        parser.add_argument("--obo", default="go-basic.obo")
        parser.add_argument("--slims", default="goslim_generic.obo")
        parser.add_argument("--gaf")
        parser.add_argument("--gene2go")
        parser.add_argument("--taxid", type=int)
        return vars(parser.parse_args(args))


# Copyright (C) 2016-present, DV Klopfenstein, H Tang, All rights reserved.
