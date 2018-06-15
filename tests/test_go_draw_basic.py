#!/usr/bin/env python
"""Test drawing GO Term diagrams through the command-line tool in go_draw.py."""

import os
import sys
from goatools.cli.gosubdag_plot import PlotCli


# pylint: disable=line-too-long

GOIDS = set(['GO:0043280', 'GO:0043281', 'GO:0043282',
             'GO:0008556', 'GO:0086037', 'GO:0005391', 'GO:0008900', 'GO:0015618'])

# STIMULUS
CMDS = [
    "--obo=i86.obo",
    "--obo=i86.obo --outfile=i86.png",
    "GO:0043280",
    "GO:0043280 GO:0043281 GO:0043282",
    "GO:0043280 GO:0043281 GO:0043282 --outfile=endopeptidase_a.png --title=Endopeptidase",
    "GO:0043280 GO:0043281 GO:0043282 -o endopeptidase_b.png -t Endopeptidase",
    "GO:0008556 GO:0086037 GO:0005391 GO:0008900 --parentcnt",
    "GO:0008556 GO:0086037 GO:0005391 GO:0008900 --parentcnt -o potas_obj_0008556.png",
    "GO:0015618 GO:0086037 GO:0005391 GO:0008900 --parentcnt -o potas_alt_0015618_d.png",
    "GO:0015618 GO:0086037 GO:0005391 GO:0008900 --childcnt --parentcnt -o potas_alt_0015618_cd.png",
    "GO:0015618 GO:0086037 GO:0005391 GO:0008900 --childcnt -o potas_alt_0015618_c.png",
]

# EXPECTED RESULTS
EXP_DICT = [
    {"obo": "i86.obo", "outfile": "go_plot.png"},
    {"obo": "i86.obo", "outfile": "i86.png"},
    {"obo":"go-basic.obo", "outfile": "go_plot.png", "GO": ["GO:0043280"]},
    {"obo":"go-basic.obo", "outfile": "go_plot.png", "GO": ["GO:0043280", "GO:0043281", "GO:0043282"]},
    {"obo":"go-basic.obo", "GO": ["GO:0043280", "GO:0043281", "GO:0043282"], "outfile": "endopeptidase_a.png", "title": "Endopeptidase"},
    {"obo":"go-basic.obo", "GO": ["GO:0043280", "GO:0043281", "GO:0043282"], "outfile": "endopeptidase_b.png", "title": "Endopeptidase"},
    {"obo":"go-basic.obo", "GO": ["GO:0008556", "GO:0086037", "GO:0005391", "GO:0008900"], "outfile": "go_plot.png"},
    {"obo":"go-basic.obo", "GO": ["GO:0008556", "GO:0086037", "GO:0005391", "GO:0008900"], "outfile": "potas_obj_0008556.png"},
    {"obo":"go-basic.obo", "GO": ["GO:0015618", "GO:0086037", "GO:0005391", "GO:0008900"], "outfile": "potas_alt_0015618_d.png"},
    {"obo":"go-basic.obo", "GO": ["GO:0015618", "GO:0086037", "GO:0005391", "GO:0008900"], "outfile": "potas_alt_0015618_cd.png"},
    {"obo":"go-basic.obo", "GO": ["GO:0015618", "GO:0086037", "GO:0005391", "GO:0008900"], "outfile": "potas_alt_0015618_c.png"},
]

EXP_SET = [
    set(),  # 0
    set(),  # 1
    set(),  # 2
    set(),  # 3
    set(),  # 4
    set(),  # 5
    set(["parentcnt"]), # 6
    set(["parentcnt"]), # 7
    set(["parentcnt"]), # 8
    set(["childcnt", "parentcnt"]), # 9
    set(["childcnt"]) # 10
]


CWD = os.getcwd()
DAT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "../../goatools/data")


def test_all(prt=sys.stdout):
    """Test all listed plotting commands."""
    # gosubdag = GoSubDag(GOIDS, get_godag())
    obj = PlotCli()
    for idx, cmdargs in enumerate(CMDS):
        prt.write("\n#######################################################\n")
        prt.write("{I}) ARGS({ARGS})\n".format(I=idx, ARGS=cmdargs))
        args = _get_arg_list(cmdargs)
        docargs = obj.get_docargs(args)
        _chk_docargs(idx, docargs)
        # prt.write("KEY                  -> VAL\n")
        # prt.write("-------------------- -- ---\n")
        # for key, val in docargs.items():
        #     prt.write("{KEY:20} -> {VAL}\n".format(KEY=key, VAL=val))
        # # assert err is None
        # # assert 'Usage' not in log

def _chk_docargs(idx, docargs):
    """Check docargs."""
    exp_d = EXP_DICT[idx]
    exp_s = EXP_SET[idx]
    print("EXP", exp_d)
    print("EXP", exp_s)
    for key_act, val_act in docargs.items():
        if key_act in set('obo'):
            val_act = os.path.basename(val_act)
        if key_act == 'outfile' and val_act[-4:] == '.png':
            val_act = os.path.basename(val_act)
        if key_act in exp_d:
            assert val_act == exp_d[key_act], "**FATAL: DICT KEY({K}) EXP({E}) ACT({V})".format(
                K=key_act, E=exp_d[key_act], V=val_act)
        elif key_act in exp_s:
            assert val_act, "**FATAL: UNKNOWN SET ELEMENT({K})".format(K=key_act)
        else:
            assert False, "UNEXPECTED"
        print("DDD KEY", key_act)

def _get_arg_list(cmdargs):
    """Turn command string into a list of args."""
    args_orig = cmdargs.split()
    args_curr = []
    for arg in args_orig:
        if arg == 'i86.obo':
            arg = os.path.join(DAT, arg)
        elif arg[-4:] == '.png':
            arg = os.path.join(CWD, arg)
        args_curr.append(arg)
    return args_orig

# class Run(object):
#     """For testing docopt of src/bin/go_plot.py."""
#
#
#
#     repo = os.path.join(os.path.dirname(os.path.abspath(__file__)), "../..")
#
#     def __init__(self):
#         self.bin = os.path.join(self.repo, "src/bin/go_plot.py")
#         self.draw_obj = PlotCli(gosubdag=self.gosubdag)
#
#     def cli(self, argstr): self.draw_obj.cli(argstr)
#     def get_outfile(self): return self.draw_obj.get_outfile()
#
#     def run_basic(self):
#         """Various plots."""
#         # One user-GO term
#         CWD = os.getcwd()
#         self.cli("GO:0043280")
#         assert os.path.exists(os.path.join(CWD, self.get_outfile()))
#         # Multiple user-GO terms
#         self.cli("GO:0043280 GO:0043281 GO:0043282")
#         assert os.path.exists(os.path.join(CWD, self.get_outfile()))
#         # Output file & Title
#         png = os.path.join(self.repo, 'endopeptidase.png')
#         docargs = {'--outfile':png}
#         self.cli("GO:0043280 GO:0043281 GO:0043282 -t Endopeptidase", **docargs)
#
#     def run_descendants(self):
#         """Annotate plots with descendants counts."""
#         # GO:0008556 is the main GO term id. GO:0015618 is an alias to GO:0008556.
#         # Print GO:0008556 when user specifies GO:0008556.
#         # Print GO:0015618 when user specifies GO:0015618, regardless that it is an alias.
#
#         # Test printing "descendantscnt" and "childcnt":
#         #
#         # descendantscnt: The total number of descendants, obtained by traversing down all
#         #     child GO term hierarchies. The count of all descendants will be printed
#         #     in the GO term box on the right side of the first text line:
#         #
#         #         GO:0086027 L09 D15 c0 d0 -> "d0" means there are 0 descendants
#         #
#         # childcnt: The number of immediate child terms for the current GO term.
#         #     Child hierarchy is not traversed. childcnt can be useful because
#         #     To prevent the GO plots from being massive and unreadable,
#         #     the default is to not plot all child terms, so knowing the total number
#         #     of immediate child terms present can give the user a better sense
#         #     of the qualities of their plot.
#         runs = [
#             "GO:0008556 GO:0086037 GO:0005391 GO:0008900 --descendantscnt -o potas_obj_0008556.png",
#             "GO:0015618 GO:0086037 GO:0005391 GO:0008900 --descendantscnt -o potas_alt_0015618_d.png",
#             "GO:0015618 GO:0086037 GO:0005391 GO:0008900 --childcnt --descendantscnt -o potas_alt_0015618_cd.png",
#             "GO:0015618 GO:0086037 GO:0005391 GO:0008900 --childcnt -o potas_alt_0015618_c.png"]
#         for arg_str in runs:
#             self.cli(arg_str)
#
#
# def test_all():
#     """Run all GO Term diagram tests."""
#     obj = Run()
#     obj.run_basic()
#     obj.run_descendants()


if __name__ == '__main__':
    test_all()
