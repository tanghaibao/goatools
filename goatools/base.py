"""Utilities used in Gene Ontology Enrichment Analyses."""

import gzip
import logging
import os
import sys
import traceback

from ftplib import FTP
from os.path import isfile

import requests

from rich.logging import RichHandler

int_types = (int,)
basestring = str


def get_logger(name: str):
    """Return a logger with a default ColoredFormatter."""
    logger = logging.getLogger(name)
    if logger.hasHandlers():
        logger.handlers.clear()
    logger.addHandler(RichHandler())
    logger.propagate = False
    logger.setLevel(logging.INFO)
    return logger


logger = get_logger("goatools")


def ungzipper(fh, blocksize=16384):
    """
    work-around to get streaming download of http://.../some.gz
    """
    import zlib

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


def download_go_basic_obo(obo="go-basic.obo", prt=sys.stdout, loading_bar=True):
    """Download Ontologies, if necessary."""
    if not isfile(obo):
        http = "http://purl.obolibrary.org/obo/go"
        if "slim" in obo:
            http = "http://www.geneontology.org/ontology/subsets"
            # http = 'http://current.geneontology.org/ontology/subsets'
        obo_remote = "{HTTP}/{OBO}".format(HTTP=http, OBO=os.path.basename(obo))
        dnld_file(obo_remote, obo, prt, loading_bar)
    else:
        if prt:
            prt.write("  EXISTS: {FILE}\n".format(FILE=obo))
    return obo


def download_ncbi_associations(gene2go="gene2go", prt=sys.stdout, loading_bar=True):
    """Download associations from NCBI, if necessary"""
    # Download: ftp://ftp.ncbi.nlm.nih.gov/gene/DATA/gene2go.gz
    gzip_file = "{GENE2GO}.gz".format(GENE2GO=gene2go)
    if not isfile(gene2go):
        file_remote = "ftp://ftp.ncbi.nlm.nih.gov/gene/DATA/{GZ}".format(
            GZ=os.path.basename(gzip_file)
        )
        dnld_file(file_remote, gene2go, prt, loading_bar)
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


def get_godag(
    fin_obo="go-basic.obo", prt=sys.stdout, loading_bar=True, optional_attrs=None
):
    """Return GODag object. Initialize, if necessary."""
    from goatools.obo_parser import GODag

    download_go_basic_obo(fin_obo, prt, loading_bar)
    return GODag(fin_obo, optional_attrs, load_obsolete=False, prt=prt)


def dnld_gaf(species_txt, prt=sys.stdout, loading_bar=True):
    """Download GAF file if necessary."""
    return dnld_gafs([species_txt], prt, loading_bar)[0]


def dnld_gafs(species_list, prt=sys.stdout, loading_bar=True):
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
        dnld_file(remove_filename, gaf_cwd, prt, loading_bar)
        fin_gafs.append(gaf_cwd)
    return fin_gafs


def http_get(url, fout=None):
    """Download a file from http. Save it in a file named by fout"""
    print("requests.get({URL}, stream=True)".format(URL=url))
    rsp = requests.get(url, stream=True)
    if rsp.status_code == 200 and fout is not None:
        with open(fout, "wb") as prt:
            for chunk in rsp:  # .iter_content(chunk_size=128):
                prt.write(chunk)
            print("  WROTE: {F}\n".format(F=fout))
    else:
        print(rsp.status_code, rsp.reason, url)
        print(rsp.content)
    return rsp


def ftp_get(fin_src, fout):
    """Download a file from an ftp server"""
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
    ftp = FTP(ftphost)  # connect to host, default port      ftp.ncbi.nlm.nih.gov
    ftp.login()  # user anonymous, passwd anonymous@
    ftp.cwd(chg_dir)  # change into "debian" directory     gene/DATA
    cmd = "RETR {F}".format(F=fin_ftp)  #                   gene2go.gz
    ftp.retrbinary(cmd, open(fout, "wb").write)  #           /usr/home/gene2go.gz
    ftp.quit()


def dnld_file(src_ftp, dst_file, prt=sys.stdout, loading_bar=True):
    """Download specified file if necessary."""
    if isfile(dst_file):
        return
    do_gunzip = src_ftp[-3:] == ".gz" and dst_file[-3:] != ".gz"
    dst_gz = "{DST}.gz".format(DST=dst_file) if do_gunzip else dst_file
    # Write to stderr, not stdout so this message will be seen when running nosetests
    cmd_msg = "get({SRC} out={DST})\n".format(SRC=src_ftp, DST=dst_gz)
    try:
        print("$ get {SRC}".format(SRC=src_ftp))
        #### wget.download(src_ftp, out=dst_gz, bar=loading_bar)
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
