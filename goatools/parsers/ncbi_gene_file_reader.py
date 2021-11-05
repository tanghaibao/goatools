"""Reads an NCBI Gene tsv file."""

from __future__ import print_function

import sys
import re
from collections import namedtuple
from collections import OrderedDict

__copyright__ = "Copyright (C) 2016-present, DV Klopfenstein, H Tang, All rights reserved."
__author__ = "DV Klopfenstein"


#pylint: disable=line-too-long,too-many-instance-attributes,unnecessary-lambda
class NCBIgeneFileReader:
    """Reads an NCBI Gene tsv file.

       Generate the NCBI gene file by following these steps:
         1) Open a browser at: https://www.ncbi.nlm.nih.gov/gene
         2) Type search text. Example:
            genetype protein coding[Properties] AND "3702"[Taxonomy ID] AND alive[property]
         3) Press the "Search" button.
         4) From the pull down menu: "Send to" -> File
    """

      # ints=None, floats=None, hdr_ex=None, log=sys.stdout):
    #def __init__(self, sep, ints, floats, hdr_ex, log):
    def __init__(self, fin, sep="\t", **kwargs_dict):
        self.log = kwargs_dict.get('log', sys.stdout)
        self.int_hdrs = [
            'tax_id', 'GeneID', 'CurrentID',  # NCBI Gene
            'start_position_on_the_genomic_accession', # NCBI Gene
            'end_position_on_the_genomic_accession',   # NCBI Gene
            'exon_count',                              # NCBI Gene
            'Start', 'start', 'End', 'end',   # Cluster
            'Len', 'len', 'Length', 'length', # cluster
            'Qty', 'qty', '# Genes']          # Cluster
        if 'ints' in kwargs_dict:
            ints = kwargs_dict['ints']
            if ints:
                self.int_hdrs.extend(ints)
            else:
                self.int_hdrs = []
        self.float_hdrs = ['Density', 'density', 'MinDensity']  # Cluster
        # These are formated for expected sorting: eg. Chr "09", "10"
        self.strpat_hdrs = {'Chr':'{:>2}', 'chromosome':'{:>2}'}
        if 'floats' in kwargs_dict:
            self.float_hdrs.extend(kwargs_dict['floats'])
        self.idxs_float = []  # run() inits proper values
        self.idxs_int = []    # run() inits proper values
        self.idxs_strpat = [] # run() inits proper values
        # Data Members used by all functions
        self.fin = fin
        self.hdr2idx = None
        self.len = 0
        self.sep = self._get_sep(fin, sep)
        self.hdr_ex = kwargs_dict.get('hdr_ex', None)
        # Data Members used by various functions
        self.ret_list = [] #             tbl2list
        self.hdrs_usr = []     # tbl2sublist tbl2list
        self.usr_max_idx = None

        # list:    Return the one item (a list of items) of interest to the user.
        # sublist: Return the items (a list of lists) of interest to the user.
        # lists:   Return all items (a list of lists) read from the tsv/csv file.
        self.fncs = {
            'list': lambda fld: self.ret_list.extend([fld[hdr_i[1]] for hdr_i in self.hdrs_usr]),
            'sublist': lambda fld: self.ret_list.append([fld[hdr_i[1]] for hdr_i in self.hdrs_usr]),
            'lists': lambda fld: self.ret_list.append(fld)
        }


    def get_h2i(self, hdrs_usr):
        """Read csv/tsv file and return specified data in a list of lists."""
        with open(self.fin) as fin_stream:
            for line in fin_stream:
                line = line.rstrip('\r\n') # chomp
                if not self.hdr2idx:
                    if self.do_hdr(line, hdrs_usr):
                        return self.hdr2idx
        return None

    def do_hdr(self, line, hdrs_usr):
        """Initialize self.h2i."""
        # If there is no header hint, consider the first line the header.
        if self.hdr_ex is None:
            self._init_hdr(line, hdrs_usr)
            return True
        # If there is a header hint, examine each beginning line until header hint is found.
        if self.hdr_ex in line:
            self._init_hdr(line, hdrs_usr)
            return True
        return False

    def run(self, fnc_name, hdrs_usr):
        """Read csv/tsv file and return specified data in a list of lists."""
        fnc = self.fncs[fnc_name]
        with open(self.fin) as fin_stream:
            for lnum, line in enumerate(fin_stream):
                line = line.rstrip('\r\n') # chomp
                # Obtain Data if headers have been collected from the first line
                if self.hdr2idx:
                    self._init_data_line(fnc, lnum, line)
                # Obtain the header
                else:
                    self.do_hdr(line, hdrs_usr)
            if self.log is not None:
                self.log.write("  {:10,} data READ:  {}\n".format(len(self.ret_list), self.fin))
        return self.ret_list, self.hdr2idx

    def get_nts(self):
        """Read csv/tsv file and return specified data in a list of lists."""
        data = []
        nt_obj = None
        with open(self.fin) as fin_stream:
            for lnum, line in enumerate(fin_stream, 1):
                try:
                    line = line.rstrip('\r\n') # chomp
                    # Obtain Data if headers have been collected from the first line
                    if nt_obj is not None:
                        flds = re.split(self.sep, line)
                        self.convert_ints_floats(flds)
                        # Aliases
                        flds[6] = self._get_list(flds[6])
                        flds[16] = [int(s) for s in self._get_list(flds[16])]
                        ntdata = nt_obj._make(flds)
                        data.append(ntdata)
                    # Obtain the header
                    else:
                        nt_obj = self._init_nt_hdr(line)
                except RuntimeError:
                    # Print headers
                    #if nt_obj is not None:
                    #  sys.stdout.write("{HDRS}\n".format(HDRS='\n'.join(nt_obj._fields)))
                    flds = re.split(self.sep, line)
                    print(len(flds), "FIELDS")
                    print(flds)
                    #raise Exception("{FIN}({LNUM}): {LINE}\n".format(
                    #  FIN=self.fin, LNUM=lnum, LINE=line))
                    # JUST SKIP LINES WITH INCOMPLETE DATA, BUT PRINT ERROR MESSAGE
                    sys.stdout.write("**ERROR: {FIN}({LNUM}): {LINE}\n".format(
                        FIN=self.fin, LNUM=lnum, LINE=line))
            if self.log is not None:
                self.log.write("  {:10,} lines READ:  {}\n".format(len(data), self.fin))
        return data

    @staticmethod
    def _get_list(valstr):
        """Return a list, given a string containing a list of values"""
        return [] if valstr == '' else [s.strip() for s in valstr.split(',')]

    def hdr_xform(self, hdrs):
        """Transform NCBI Gene header fields into valid namedtuple fields."""
        xform = []
        hdrs = self.replace_nulls(hdrs)
        for hdr in hdrs:
            hdr = hdr.replace('.', '_')
            hdr = hdr.replace(' ', '_')
            hdr = hdr.replace('#', 'N')
            hdr = hdr.replace('-', '_')
            hdr = hdr.replace('"', '')
            xform.append(hdr)
        return xform

    def _init_nt_hdr(self, line):
        """Convert headers into valid namedtuple fields."""
        line = line.replace('.', '_')
        line = line.replace(' ', '_')
        line = line.replace('#', 'N')
        line = line.replace('-', '_')
        line = line.replace('"', '')
        #line = re.sub(r"_$", r"", line)
        hdrs = re.split(self.sep, line)
        if '' in hdrs:
            hdrs = NCBIgeneFileReader.replace_nulls(hdrs)
        # Init indexes which will be converted to int or float
        self.idxs_int = [idx for idx, hdr in enumerate(hdrs) if hdr in self.int_hdrs]
        self.idxs_float = [idx for idx, hdr in enumerate(hdrs) if hdr in self.float_hdrs]
        if hdrs[6] != 'Aliases':
            raise RuntimeError('**FATAL: BAD HEADER LINE: {LINE}'.format(LINE=line))
        return namedtuple('ntncbi', ' '.join(hdrs))

    @staticmethod
    def _get_sep(fin, sep):
        """Uses extension(.tsv, .csv) to determine separator."""
        if '.tsv' in fin:
            return r'\t'
        if '.csv' in fin:
            return r','
        return sep

    @staticmethod
    def replace_nulls(hdrs):
        """Replace '' in hdrs."""
        ret = []
        idx = 0
        for hdr in hdrs:
            if hdr == '':
                ret.append("no_hdr{}".format(idx))
            else:
                ret.append(hdr)
        return ret

    def _init_data_line(self, fnc, lnum, line):
        """Process Data line."""
        fld = re.split(self.sep, line)
        # Lines may contain different numbers of items.
        # The line should have all columns requested by the user.
        if self.usr_max_idx < len(fld):
            self.convert_ints_floats(fld)
            fnc(fld)
        else:
            for fld in enumerate(zip(self.hdr2idx.keys(), fld)):
                print(fld)
            for hdr in self.hdrs_usr:
                print(hdr)
            print('# ITEMS ON A LINE:', len(fld))
            print('MAX USR IDX:', self.usr_max_idx)
            raise Exception("ERROR ON LINE {} IN {}".format(lnum+1, self.fin))

    def convert_ints_floats(self, flds):
        """Convert strings to ints and floats, if so specified."""
        for idx in self.idxs_float:
            flds[idx] = float(flds[idx])
        for idx in self.idxs_int:
            dig = flds[idx]
            #print 'idx={} ({}) {}'.format(idx, flds[idx], flds) # DVK
            flds[idx] = int(flds[idx]) if dig.isdigit() else dig
        for idx in self.idxs_strpat:
            hdr = self.hdr2idx.items()[idx][0]
            pat = self.strpat_hdrs[hdr]
            flds[idx] = pat.format(flds[idx])

    def _init_hdr(self, line, hdrs_usr):
        """Initialize self.hdr2idx, self.len, self.idxs_float, and self.idxs_int"""
        self.hdr2idx = OrderedDict([(v.strip(), i) for i, v in enumerate(re.split(self.sep, line))])
        self.len = len(self.hdr2idx)
        # If user is requesting specific data fields...
        if hdrs_usr is not None:
            # Loop through the user headers
            for usr_hdr in hdrs_usr:
                # If the user header is contained in the file....
                if usr_hdr in self.hdr2idx:
                    # Add the user header and the field index to a list
                    self.hdrs_usr.append([usr_hdr, self.hdr2idx[usr_hdr]])
                else:
                    raise Exception("NO COLUMN({}) FOUND:\n  HDR={}\n".format(
                        hdrs_usr, '\n  HDR='.join(self.hdr2idx.keys())))
        usr_hdrs = [E[0] for E in self.hdrs_usr] if self.hdrs_usr else self.hdr2idx
        self._init_idxs_float(usr_hdrs)
        self._init_idxs_int(usr_hdrs)
        self._init_idxs_strpat(usr_hdrs)
        self.usr_max_idx = max(E[1] for E in self.hdrs_usr) if self.hdrs_usr else len(self.hdr2idx)-1

    def _init_idxs_float(self, usr_hdrs):
        """List of indexes whose values will be floats."""
        self.idxs_float = [
            Idx for Hdr, Idx in self.hdr2idx.items() if Hdr in usr_hdrs and Hdr in self.float_hdrs]

    def _init_idxs_int(self, usr_hdrs):
        """List of indexes whose values will be ints."""
        self.idxs_int = [
            Idx for Hdr, Idx in self.hdr2idx.items() if Hdr in usr_hdrs and Hdr in self.int_hdrs]

    def _init_idxs_strpat(self, usr_hdrs):
        """List of indexes whose values will be strings."""
        strpat = self.strpat_hdrs.keys()
        self.idxs_strpat = [
            Idx for Hdr, Idx in self.hdr2idx.items() if Hdr in usr_hdrs and Hdr in strpat]


  # Copyright (C) 2016-present, DV Klopfenstein, H Tang, All rights reserved.
