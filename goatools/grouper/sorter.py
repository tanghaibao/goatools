"""Sorts GO IDs or user-provided sections containing GO IDs."""

__copyright__ = "Copyright (C) 2016-2019, DV Klopfenstein, H Tang, All rights reserved."
__author__ = "DV Klopfenstein"

import sys
import collections as cx
from goatools.wr_tbl import prt_txt
from goatools.grouper.sorter_nts import SorterNts
from goatools.grouper.sorter_gos import SorterGoIds
from goatools.grouper.wr_sections import WrSectionsTxt


class Sorter(object):
    """Sorts GO IDs or user-provided sections containing GO IDs.

       User GO IDs grouped under header GO IDs are not sorted by the Grouper class.
       Sort both user GO IDs in a group and header GO IDs across groups with these:

       S: use_sections
       s: section_sortby (T=True, F=False, S=lambda sort function)
       h: hdrgo_sortby Sorts hdr GO IDs
       u: sortby       Sorts user GO IDs
       P: hdrgo_prt    If True, Removes GO IDs used as GO group headers; Leaves list in
                       sorted order, but removes header GO IDs which are not user GO IDs.

               rm_h     hdr_sort      usr_sort       S  s  h  u  p
                ---     ------------  ------------  -- -- -- -- --
       case 1:   NO     hdrgo_sortby  usrgo_sortby   N  T  H  U  T
       case 2:  YES     hdrgo_sortby  usrgo_sortby   N  T  H  U  F
       case 3:   NO     section_order usrgo_sortby   S  F  -  U  T
       case 4:  YES     section_order usrgo_sortby   S  F  -  U  F
       case 5:  YES     |<--- section_sortby --->|   S  S  -  -  -

                                        |print|
         sec usesec prthdr prtsec 1d 2d hdr usr
        ---- ------ ------ ------ -- -- --- ---
        none      -   true      -  y  . hdr usr A
        none      -  false      -  y  . ... usr B

         yes  False   true      -  y  . hdr usr A
         yes  False  false      -  y  . ... usr B

         yes   True   True  False  .  y hdr usr
         yes   True  False  False  .  y ... usr
    """

    # Keywords for creating desc2nts
    keys_nts = set(["hdrgo_prt", "section_prt", "top_n", "use_sections"])

    def __init__(self, grprobj, **kws):
        # Keyword arguments:
        _sortby = kws.get('sortby')
        _hdrgo_sortby = kws.get('hdrgo_sortby')
        _section_sortby = kws.get('section_sortby')
        # GO IDs are grouped, but not yet sorted
        # print('SSSSSSSSSSS Sorter(sortby={} hdrgo_sortby={}, section_sortby={}'.format(
        #     _sortby, _hdrgo_sortby, _section_sortby))
        self.grprobj = grprobj
        # SorterGoIds can return either a 2-D list of sorted GO IDs or a flat sorted GO list
        self.sortgos = SorterGoIds(grprobj, _sortby, _hdrgo_sortby)
        self.sectobj = SorterNts(self.sortgos, _section_sortby) if grprobj.hdrobj.sections else None

    def prt_gos(self, prt=sys.stdout, **kws_usr):
        """Sort user GO ids, grouped under broader GO terms or sections. Print to screen."""
        # deprecated
        # Keyword arguments (control content): hdrgo_prt section_prt use_sections
        # desc2nts contains: (sections hdrgo_prt sortobj) or (flat hdrgo_prt sortobj)
        desc2nts = self.get_desc2nts(**kws_usr)
        # Keyword arguments (control print format): prt prtfmt
        self.prt_nts(desc2nts, prt, kws_usr.get('prtfmt'))
        return desc2nts

    def prt_nts(self, desc2nts, prt=sys.stdout, prtfmt=None):
        """Print grouped and sorted GO IDs."""
        # deprecated
        # Set print format string
        if prtfmt is None:
            prtfmt = "{{hdr1usr01:2}} {FMT}\n".format(FMT=self.grprobj.gosubdag.prt_attr['fmt'])
        # 1-D: data to print is a flat list of namedtuples
        if 'flat' in desc2nts:
            prt_txt(prt, desc2nts['flat'], prtfmt=prtfmt)
        # 2-D: data to print is a list of [(section, nts), ...
        else:
            WrSectionsTxt.prt_sections(prt, desc2nts['sections'], prtfmt)

    def get_desc2nts(self, **kws_usr):
        """Return grouped, sorted namedtuples in either format: flat, sections."""
        # desc2nts contains: (sections hdrgo_prt sortobj) or (flat hdrgo_prt sortobj)
        # keys_nts: hdrgo_prt section_prt top_n use_sections
        kws_nts = {k:v for k, v in kws_usr.items() if k in self.keys_nts}
        return self.get_desc2nts_fnc(**kws_nts)

    def get_desc2nts_fnc(self, hdrgo_prt=True, section_prt=None,
                         top_n=None, use_sections=True):
        """Return grouped, sorted namedtuples in either format: flat, sections."""
        # RETURN: flat list of namedtuples
        nts_flat = self.get_nts_flat(hdrgo_prt, use_sections)
        if nts_flat:
            flds = nts_flat[0]._fields
            if not use_sections:
                return {'sortobj':self, 'flat' : nts_flat, 'hdrgo_prt':hdrgo_prt, 'flds':flds,
                        'num_items':len(nts_flat), 'num_sections':1}
            else:
                return {'sortobj':self,
                        'sections' : [(self.grprobj.hdrobj.secdflt, nts_flat)],
                        'hdrgo_prt':hdrgo_prt,
                        'flds':flds,
                        'num_items':len(nts_flat), 'num_sections':1}
        # print('FFFF Sorter:get_desc2nts_fnc: nts_flat is None')
        # RETURN: 2-D list [(section_name0, namedtuples0), (section_name1, namedtuples1), ...
        #     kws: top_n hdrgo_prt section_sortby
        # Over-ride hdrgo_prt depending on top_n value
        assert top_n is not True and top_n is not False, \
            "top_n({T}) MUST BE None OR AN int".format(T=top_n)
        assert self.sectobj is not None, "SECTIONS OBJECT DOES NOT EXIST"
        sec_sb = self.sectobj.section_sortby
        # Override hdrgo_prt, if sorting by sections or returning a subset of GO IDs in section
        hdrgo_prt_curr = hdrgo_prt is True
        if sec_sb is True or (sec_sb is not False and sec_sb is not None) or top_n is not None:
            hdrgo_prt_curr = False
        # print('GGGG Sorter:get_desc2nts_fnc: hdrgo_prt_curr({}) sec_sb({}) top_n({})'.format(
        #     hdrgo_prt_curr, sec_sb, top_n))
        nts_section = self.sectobj.get_sorted_nts_keep_section(hdrgo_prt_curr)
        # print('HHHH Sorter:get_desc2nts_fnc: nts_section')
        # Take top_n in each section, if requested
        if top_n is not None:
            nts_section = [(s, nts[:top_n]) for s, nts in nts_section]
            if section_prt is None:
                nts_flat = self.get_sections_flattened(nts_section)
                flds = nts_flat[0]._fields if nts_flat else []
                return {'sortobj':self, 'flat' : nts_flat, 'hdrgo_prt':hdrgo_prt_curr, 'flds':flds,
                        'num_items':len(nts_flat), 'num_sections':1}
        # Send flat list of sections nts back, as requested
        if section_prt is False:
            nts_flat = self.get_sections_flattened(nts_section)
            flds = nts_flat[0]._fields if nts_flat else []
            return {'sortobj':self, 'flat' : nts_flat, 'hdrgo_prt':hdrgo_prt_curr, 'flds':flds,
                    'num_items':len(nts_flat),
                    'num_sections':len(nts_section)}
        # Send 2-D sections nts back
        # print('IIII Sorter:get_desc2nts_fnc: nts_section')
        flds = nts_section[0][1][0]._fields if nts_section else []
        return {'sortobj':self, 'sections' : nts_section, 'hdrgo_prt':hdrgo_prt_curr, 'flds':flds,
                'num_items':sum(len(nts) for _, nts in nts_section),
                'num_sections':len(nts_section)}

    @staticmethod
    def get_sections_flattened(section_nts):
        """Convert [(section0, nts0), (section1, nts1), ... to [*nts0, *nts1, ..."""
        nt_flds = list(section_nts[0][1][0]._fields)
        # Flatten section_nts 2-D list
        if 'section' in nt_flds:
            return [nt for _, nts in section_nts for nt in nts]
        # Flatten section_nts 2-D list, and add sections to each namedtuple
        nt_flds.append('section')
        nts_flat = []
        ntobj = cx.namedtuple("Nt", " ".join(nt_flds))
        for section_name, nts in section_nts:
            for nt_go in nts:
                vals = list(nt_go) + [section_name]
                nts_flat.append(ntobj._make(vals))
        return nts_flat

    def get_nts_flat(self, hdrgo_prt=True, use_sections=True):
        """Return a flat list of sorted nts."""
        # Either there are no sections OR we are not using them
        if self.sectobj is None or not use_sections:
            return self.sortgos.get_nts_sorted(
                hdrgo_prt,
                hdrgos=self.grprobj.get_hdrgos(),
                hdrgo_sort=True)
        if not use_sections:
            return self.sectobj.get_sorted_nts_omit_section(hdrgo_prt, hdrgo_sort=True)
        return None

    @staticmethod
    def get_fields(desc2nts):
        """Return grouped, sorted namedtuples in either format: flat, sections."""
        if 'flat' in desc2nts:
            nts_flat = desc2nts.get('flat')
            if nts_flat:
                return nts_flat[0]._fields
        if 'sections' in desc2nts:
            nts_sections = desc2nts.get('sections')
            if nts_sections:
                return nts_sections[0][1][0]._fields


# Copyright (C) 2016-2019, DV Klopfenstein, H Tang, All rights reserved.
