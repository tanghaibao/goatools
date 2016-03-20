# Stolen from brentp:
# <https://github.com/brentp/toolshed/blob/master/toolshed/files.py>

import bz2
import gzip
import sys
import urllib


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
