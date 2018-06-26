"""For getting colors based on grouping characteristics."""

__copyright__ = "Copyright (C) 2016-2018, DV Klopfenstein, H Tang, All rights reserved."
__author__ = "DV Klopfenstein"


class GrouperColors(object):
    """Groups the user GO ids under other GO IDs acting as headers for the GO groups."""

    # hdrcol = '#1f0954'  # dark indigo
    hdrcol_all = '#152eff'  # vivid blue

    def __init__(self, grprobj):
        self.grprobj = grprobj
        self.usrgos = grprobj.usrgos
        self.hdrgos_actual = grprobj.get_hdrgos()

    def get_bordercolor(self):
        """Get bordercolor based on hdrgos and usergos."""
        hdrgos_all = self.grprobj.hdrobj.get_hdrgos()
        hdrgos_unused = hdrgos_all.difference(self.hdrgos_actual)
        go2bordercolor = {}
        # hdrgos that went unused
        for hdrgo in hdrgos_unused:
            go2bordercolor[hdrgo] = self.hdrcol_all
        # hdrgos used in this grouping that are NOT usrgos
        for hdrgo in self.grprobj.hdrgo2usrgos.keys():
            go2bordercolor[hdrgo] = self.hdrcol_all
        # hdrgos used in this grouping that ARE usrgos
        for hdrgo in self.grprobj.hdrgo_is_usrgo:
            go2bordercolor[hdrgo] = 'blue'
        # usrgos which are NOT hdrgos
        usrgos_rem = self.grprobj.usrgos.difference(self.grprobj.hdrgo_is_usrgo)
        for usrgo in usrgos_rem:
            go2bordercolor[usrgo] = '#029386'  # teal
        # print("{N:5} hdrgos actual".format(N=len(self.hdrgos_actual)))
        # print("{N:5} hdrgos unused".format(N=len(hdrgos_unused)))
        # print("{N:5} hdrgos only       BLACK".format(N=len(self.grprobj.hdrgo2usrgos.keys())))
        # print("{N:5} usrgos".format(N=len(self.grprobj.usrgos)))
        # print("{N:5} usrgos AND hdrgos BLUE".format(N=len(self.grprobj.hdrgo_is_usrgo)))
        # print("{N:5} usrgos Only".format(N=len(usrgos_rem)))
        return go2bordercolor

    def get_go2color_users(self,
                           usrgo_color='#feffa3', # yellow
                           hdrusrgo_color='#d4ffea', # green
                           hdrgo_color='#eee6f6'): # purple
        """Get go2color for GO DAG plots."""
        go2color = {}
        # Color user GO IDs
        for goid in self.usrgos:
            go2color[goid] = usrgo_color
        # Color header GO IDs. Headers which are also GO IDs get their own color.
        for goid_hdr in self.hdrgos_actual:
            go2color[goid_hdr] = hdrusrgo_color if goid_hdr in self.usrgos else hdrgo_color
        return go2color

# Copyright (C) 2016-2018, DV Klopfenstein, H Tang, All rights reserved.
