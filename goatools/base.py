"""Utilities used in Gene Ontology Enrichment Analyses."""
# Stolen from brentp:
# <https://github.com/brentp/toolshed/blob/master/toolshed/files.py>

import os
import os.path as op
import sys
import bz2
import gzip
import urllib
import wget


if sys.version_info[0] < 3:
    int_types = (int, long)
    urlopen = urllib.urlopen
else:
    int_types = (int,)
    basestring = str
    from urllib.request import urlopen


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
    if isinstance(f, int_types):
        return nopen(sys.argv[f], mode)

    if not isinstance(f, basestring):
        return f
    if f.startswith("|"):
        # using shell explicitly makes things like process substitution work:
        # http://stackoverflow.com/questions/7407667/python-subprocess-subshells-and-redirection
        # use sys.stderr so we dont have to worry about checking it...
        p = Popen(f[1:], stdout=PIPE, stdin=PIPE,
                  stderr=sys.stderr if mode == "r" else PIPE,
                  shell=True, bufsize=-1, # use system default for buffering
                  preexec_fn=prefunc,
                  close_fds=False, executable=os.environ.get('SHELL'))
        if sys.version_info[0] > 2:
            import io
            p.stdout = io.TextIOWrapper(p.stdout)
            p.stdin = io.TextIOWrapper(p.stdin)
            if mode != "r":
                p.stderr = io.TextIOWrapper(p.stderr)

        if mode and mode[0] == "r":
            return process_iter(p, f[1:])
        return p

    if f.startswith(("http://", "https://", "ftp://")):
        fh = urlopen(f)
        if f.endswith(".gz"):
            return ungzipper(fh)
        if sys.version_info[0] < 3:
            return fh
        import io
        return io.TextIOWrapper(fh)
    f = op.expanduser(op.expandvars(f))
    if f.endswith((".gz", ".Z", ".z")):
        fh = gzip.open(f, mode)
        if sys.version_info[0] < 3:
            return fh
        import io
        return io.TextIOWrapper(fh)
    elif f.endswith((".bz", ".bz2", ".bzip2")):
        fh = bz2.BZ2File(f, mode)
        if sys.version_info[0] < 3:
            return fh
        import io
        return io.TextIOWrapper(fh)

    return {"r": sys.stdin, "w": sys.stdout}[mode[0]] if f == "-" \
         else open(f, mode)


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
    if not os.path.isfile(obo):
        http = "http://purl.obolibrary.org/obo/go"
        if "slim" in obo:
            http = "http://www.geneontology.org/ontology/subsets"
        obo_remote = "{HTTP}/{OBO}".format(HTTP=http, OBO=os.path.basename(obo))
        dnld_file(obo_remote, obo, prt, loading_bar)
    else:
        if prt is not None:
            prt.write("  EXISTS: {FILE}\n".format(FILE=obo))
    return obo

def download_ncbi_associations(gene2go="gene2go", prt=sys.stdout, loading_bar=True):
    """Download associations from NCBI, if necessary"""
    # Download: ftp://ftp.ncbi.nlm.nih.gov/gene/DATA/gene2go.gz
    gzip_file = "{GENE2GO}.gz".format(GENE2GO=gene2go)
    if not os.path.isfile(gene2go):
        file_remote = "ftp://ftp.ncbi.nlm.nih.gov/gene/DATA/{GZ}".format(
            GZ=os.path.basename(gzip_file))
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

def get_godag(fin_obo="go-basic.obo", prt=sys.stdout, loading_bar=True):
    """Return GODag object. Initialize, if necessary."""
    from goatools.obo_parser import GODag
    download_go_basic_obo(fin_obo, prt, loading_bar)
    return GODag(fin_obo)

def get_gaf_name(species):
    """Given a species (eg goa_human, mgi, fb), return filename of GAF file."""
    gaf_pats = {
        'gas':"gene_association.{S}",
        'goa':"{S}.gaf"}
    # Example species text: goa_human mgi fb
    gaf_key = 'goa' if species[:4] == "goa_" else 'gas'
    # Return Examples: goa_human.gaf gene_association.mgi gene_association.fb
    return gaf_pats[gaf_key].format(S=species)

def dnld_gaf(species_txt, prt=sys.stdout, loading_bar=True):
    """Download GAF file if necessary."""
    return dnld_gafs([species_txt], prt, loading_bar)[0]

def dnld_gafs(species_list, prt=sys.stdout, loading_bar=True):
    """Download GAF files if necessary."""
    # Example GAF files:
    #   http://geneontology.org/gene-associations/gene_association.mgi.gz
    #   http://geneontology.org/gene-associations/gene_association.fb.gz
    #   http://geneontology.org/gene-associations/goa_human.gaf.gz
    #   NA: http://geneontology.org/gene-associations/gene_association.goa_human.gz
    http = "http://geneontology.org/gene-associations"
    # There are two filename patterns for gene associations on geneontology.org
    fin_gafs = []
    cwd = os.getcwd()
    for species_txt in species_list: # e.g., goa_human mgi fb
        gaf_base = get_gaf_name(species_txt) # goa_human.gaf
        gaf_cwd = os.path.join(cwd, gaf_base) # {CWD}/goa_human.gaf
        wget_cmd = "{HTTP}/{GAF}.gz".format(HTTP=http, GAF=gaf_base)
        dnld_file(wget_cmd, gaf_cwd, prt, loading_bar)
        fin_gafs.append(gaf_cwd)
    return fin_gafs

def dnld_file(src_ftp, dst_file, prt=sys.stdout, loading_bar=True):
    """Download specified file if necessary."""
    if os.path.isfile(dst_file):
        return
    do_gunzip = src_ftp[-3:] == '.gz' and dst_file[-3:] != '.gz'
    dst_wget = "{DST}.gz".format(DST=dst_file) if do_gunzip else dst_file
    # Write to stderr, not stdout so this message will be seen when running nosetests
    wget_msg = "wget.download({SRC} out={DST})\n".format(SRC=src_ftp, DST=dst_wget)
    sys.stderr.write("  {WGET}".format(WGET=wget_msg))
    if loading_bar:
        loading_bar = wget.bar_adaptive
    try:
        wget.download(src_ftp, out=dst_wget, bar=loading_bar)
        if do_gunzip:
            if prt is not None:
                prt.write("  gunzip {FILE}\n".format(FILE=dst_wget))
            gzip_open_to(dst_wget, dst_file)
    except IOError as errmsg:
        import traceback
        traceback.print_exc()
        sys.stderr.write("**FATAL cmd: {WGET}".format(WGET=wget_msg))
        sys.stderr.write("**FATAL msg: {ERR}".format(ERR=str(errmsg)))
        sys.exit(1)

def gzip_open_to(fin_gz, fout):
    """Unzip a file.gz file."""
    with gzip.open(fin_gz, 'rb') as zstrm:
        with  open(fout, 'wb') as ostrm:
            ostrm.write(zstrm.read())
    assert os.path.isfile(fout), "COULD NOT GUNZIP({G}) TO FILE({F})".format(G=fin_gz, F=fout)
    os.remove(fin_gz)

# Copyright (C) 2013-2018, B Pedersen, et al. All rights reserved."
