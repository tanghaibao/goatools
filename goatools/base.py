"""Utilities used in Gene Ontology Enrichment Analyses."""

import bz2
import gzip
import io
import logging
import os
import os.path as op
import sys
import traceback
import zlib

from os.path import isfile
from subprocess import PIPE, Popen
from urllib.request import urlopen

import requests

from ftpretty import ftpretty
from rich.logging import RichHandler


def get_logger(name: str):
    """
    Return a logger with a default ColoredFormatter.
    """
    log = logging.getLogger(name)
    if log.hasHandlers():
        log.handlers.clear()
    log.addHandler(RichHandler())
    log.propagate = False
    log.setLevel(logging.INFO)
    return log


logger = get_logger("goatools")


def nopen(f, mode="r"):
    r"""
    open a file that's gzipped or return stdin for '-'
    if f is a number, the result of nopen(sys.argv[f]) is returned.
    >>> nopen('-') == sys.stdin, nopen('-', 'w') == sys.stdout
    (True, True)
    >>> nopen(sys.argv[0])
    <...file...>
    # expands user and vars ($HOME)
    >>> nopen("~/.bashrc").name == nopen("$HOME/.bashrc").name
    True
    # an already open file.
    >>> nopen(open(sys.argv[0]))
    <...file...>
    >>> nopen(0)
    <...file...>
    Or provide nicer access to Popen.stdout
    >>> files = list(nopen("|ls"))
    >>> assert 'setup.py\n' in files or b'setup.py\n' in files, files
    """
    if isinstance(f, int):
        return nopen(sys.argv[f], mode)

    if not isinstance(f, str):
        return f
    if f.startswith("|"):
        # using shell explicitly makes things like process substitution work:
        # http://stackoverflow.com/questions/7407667/python-subprocess-subshells-and-redirection
        # use sys.stderr so we dont have to worry about checking it...
        p = Popen(
            f[1:],
            stdout=PIPE,
            stdin=PIPE,
            stderr=sys.stderr if mode == "r" else PIPE,
            shell=True,
            bufsize=-1,  # use system default for buffering
            close_fds=False,
            executable=os.environ.get("SHELL"),
        )
        p.stdout = io.TextIOWrapper(p.stdout)
        p.stdin = io.TextIOWrapper(p.stdin)
        if mode != "r":
            p.stderr = io.TextIOWrapper(p.stderr)

        return p

    if f.startswith(("http://", "https://", "ftp://")):
        fh = urlopen(f)
        if f.endswith(".gz"):
            return ungzipper(fh)
        return io.TextIOWrapper(fh)
    f = op.expanduser(op.expandvars(f))
    if f.endswith((".gz", ".Z", ".z")):
        fh = gzip.open(f, mode)
        return io.TextIOWrapper(fh)
    elif f.endswith((".bz", ".bz2", ".bzip2")):
        fh = bz2.BZ2File(f, mode)
        return io.TextIOWrapper(fh)

    return (
        {"r": sys.stdin, "w": sys.stdout}[mode[0]]
        if f == "-"
        else open(f, mode, encoding="utf-8")
    )


def ungzipper(fh, blocksize=16384):
    """
    work-around to get streaming download of http://.../some.gz
    """
    uzip = zlib.decompressobj(16 + zlib.MAX_WBITS)
    data = uzip.decompress(fh.read(blocksize)).split("\n")

    while len(data[0]):
        # last chunk might not be a full line.
        save = data.pop()
        for line in data:
            yield line
        data = uzip.decompress(fh.read(blocksize)).split("\n")
        # first line is prepended with saved chunk from end of last set.
        data[0] = save + data[0]


def download_go_basic_obo(obo="go-basic.obo", prt=sys.stdout):
    """Download Ontologies, if necessary."""
    if not isfile(obo):
        http = "http://purl.obolibrary.org/obo/go"
        if "slim" in obo:
            http = "http://www.geneontology.org/ontology/subsets"
        obo_remote = f"{http}/{op.basename(obo)}"
        dnld_file(obo_remote, obo, prt)
    else:
        if prt:
            prt.write("  EXISTS: {FILE}\n".format(FILE=obo))
    return obo


def download_ncbi_associations(gene2go="gene2go", prt=sys.stdout):
    """Download associations from NCBI, if necessary"""
    # Download: ftp://ftp.ncbi.nlm.nih.gov/gene/DATA/gene2go.gz
    gzip_file = "{GENE2GO}.gz".format(GENE2GO=gene2go)
    if not isfile(gene2go):
        file_remote = f"ftp://ftp.ncbi.nlm.nih.gov/gene/DATA/{op.basename(gzip_file)}"
        dnld_file(file_remote, gene2go, prt)
    else:
        if prt is not None:
            prt.write("  EXISTS: {FILE}\n".format(FILE=gene2go))
    return gene2go


def gunzip(gzip_file, file_gunzip=None):
    """Unzip .gz file. Return filename of unzipped file."""
    if file_gunzip is None:
        file_gunzip = os.path.splitext(gzip_file)[0]
        gzip_open_to(gzip_file, file_gunzip)
        return file_gunzip


def get_godag(fin_obo="go-basic.obo", prt=sys.stdout, optional_attrs=None):
    """Return GODag object. Initialize, if necessary."""
    from .obo_parser import GODag

    download_go_basic_obo(fin_obo, prt)
    return GODag(fin_obo, optional_attrs, load_obsolete=False, prt=prt)


def dnld_gaf(species_txt, prt=sys.stdout):
    """Download GAF file if necessary."""
    return dnld_gafs([species_txt], prt)[0]


def dnld_gafs(species_list, prt=sys.stdout):
    """Download GAF files if necessary."""
    # Example GAF files in  http://current.geneontology.org/annotations/:
    #   http://current.geneontology.org/annotations/mgi.gaf.gz
    #   http://current.geneontology.org/annotations/fb.gaf.gz
    #   http://current.geneontology.org/annotations/goa_human.gaf.gz
    http = "http://current.geneontology.org/annotations"
    # There are two filename patterns for gene associations on geneontology.org
    fin_gafs = []
    cwd = os.getcwd()
    for species_txt in species_list:  # e.g., goa_human mgi fb
        gaf_base = "{ABC}.gaf".format(ABC=species_txt)  # goa_human.gaf
        gaf_cwd = os.path.join(cwd, gaf_base)  # {CWD}/goa_human.gaf
        remove_filename = "{HTTP}/{GAF}.gz".format(HTTP=http, GAF=gaf_base)
        dnld_file(remove_filename, gaf_cwd, prt)
        fin_gafs.append(gaf_cwd)
    return fin_gafs


def http_get(url, fout=None):
    """Download a file from http. Save it in a file named by fout"""
    print("requests.get({URL}, stream=True)".format(URL=url))
    rsp = requests.get(url, stream=True, timeout=10)
    if rsp.status_code == 200 and fout is not None:
        with open(fout, "wb") as prt:
            for chunk in rsp:  # .iter_content(chunk_size=128):
                prt.write(chunk)
            print("  WROTE: {F}\n".format(F=fout))
    else:
        print(rsp.status_code, rsp.reason, url)
        print(rsp.content)
    return rsp


def ftp_get(fin_src: str, fout: str):
    """
    Download a file from an ftp server, e.g., ftp://ftp.ncbi.nlm.nih.gov/gene/DATA/gene2go.gz
    """
    assert fin_src[:6] == "ftp://", fin_src
    dir_full, fin_ftp = os.path.split(fin_src[6:])
    pt0 = dir_full.find("/")
    assert pt0 != -1, pt0
    ftphost = dir_full[:pt0]
    chg_dir = dir_full[pt0 + 1 :]
    print(
        "FTP RETR {HOST} {DIR} {SRC} -> {DST}".format(
            HOST=ftphost, DIR=chg_dir, SRC=fin_ftp, DST=fout
        )
    )
    ftp = ftpretty(ftphost, "anonymous", "anonymous@")
    ftp.get(chg_dir + "/" + fin_ftp, fout)


def dnld_file(src_ftp, dst_file, prt=sys.stdout):
    """Download specified file if necessary."""
    if isfile(dst_file):
        return
    do_gunzip = src_ftp[-3:] == ".gz" and dst_file[-3:] != ".gz"
    dst_gz = "{DST}.gz".format(DST=dst_file) if do_gunzip else dst_file
    # Write to stderr, not stdout so this message will be seen when running nosetests
    cmd_msg = "get({SRC} out={DST})\n".format(SRC=src_ftp, DST=dst_gz)
    try:
        print("$ get {SRC}".format(SRC=src_ftp))
        if src_ftp[:4] == "http":
            http_get(src_ftp, dst_gz)
        else:
            ftp_get(src_ftp, dst_gz)
        if do_gunzip:
            if prt is not None:
                prt.write("$ gunzip {FILE}\n".format(FILE=dst_gz))
            gzip_open_to(dst_gz, dst_file)
    except IOError as errmsg:
        traceback.print_exc()
        logger.fatal("cmd: %s", cmd_msg)
        logger.fatal("msg: %s", str(errmsg))
        sys.exit(1)


def gzip_open_to(fin_gz, fout):
    """Unzip a file.gz file."""
    try:
        with gzip.open(fin_gz, "rb") as zstrm:
            with open(fout, "wb") as ostrm:
                ostrm.write(zstrm.read())
    # pylint: disable=broad-except
    except Exception as errmsg:
        logger.error("COULD NOT GUNZIP(%s) TO FILE(%s)", fin_gz, fout)
        traceback.print_exc()
        logger.fatal("msg: %s", str(errmsg))
        sys.exit(1)
    os.remove(fin_gz)


# Copyright (C) 2013-present, B Pedersen, et al. All rights reserved."
