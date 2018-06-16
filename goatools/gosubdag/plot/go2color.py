"""Manages GO Term fill colors and bordercolors."""

__copyright__ = "Copyright (C) 2016-2018, DV Klopfenstein, H Tang, All rights reserved."
__author__ = "DV Klopfenstein"

import collections as cx
from goatools.gosubdag.utils import get_kwargs


class Go2Color(object):
    """Manages GO Term fill colors and bordercolors."""

    kws_dct = set(['go2color', 'go2bordercolor', 'dflt_bordercolor', 'key2col'])

    alpha2col = cx.OrderedDict([
        # GOEA GO terms that are significant
        (0.005, 'mistyrose'),
        (0.010, 'moccasin'),
        (0.050, 'lemonchiffon1'),
        # GOEA GO terms that are not significant
        (1.000, 'grey95'),
    ])

    key2col = {
        'level_01': '#f1fbfd', # Very light blue
        'go_sources': '#ffffe4', # xkcd off white
        'go_leafchild': '#eae8e8', # 'Grey palette' http://www.color-hex.com/color-palette/45491
    }


    def __init__(self, gosubdag, objgoea=None, **kws):
        # kws: go2color go2bordercolor dflt_bordercolor
        self.kws = get_kwargs(kws, self.kws_dct, None)
        # Use default key coloring if user does not specify to turn it off (key2col=False)
        if 'key2col' not in kws or kws['key2col']:
            self.kws['key2col'] = self.key2col
        self.gosubdag = gosubdag  # GoSubDag
        self.objgoea = objgoea    # GoeaResults
        self.dflt_bordercolor = kws.get('dflt_bordercolor', 'mediumseagreen')
        self.go2color = self.init_goid2color()
        self.go2bordercolor = kws.get('go2bordercolor', {})
        self._init_equiv()

    def get_bordercolor(self, goid):
        """Return GO Term border color."""
        return self.go2bordercolor.get(goid, self.dflt_bordercolor)

    def init_goid2color(self):
        """Set colors of GO terms."""
        goid2color = {}
        # 1. User-specified colors for each GO term
        if 'go2color' in self.kws:
            for goid, color in self.kws['go2color'].items():
                goid2color[goid] = color
        # 2. colors based on p-value override colors based on source GO
        if self.objgoea is not None:
            self.objgoea.set_goid2color_pval(goid2color)
        key2color = self.kws.get('key2col')
        if key2color is not None:
            # 3. Default GO source color
            if 'go_sources' in key2color:
                color = key2color['go_sources']
                go2obj = self.gosubdag.go2obj
                for goid in self.gosubdag.go_sources:
                    if goid not in goid2color:
                        goobj = go2obj[goid]
                        goid2color[goobj.id] = color
                        goid2color[goid] = color
            # 4. Level-01 GO color
            if 'level_01' in key2color:
                color = key2color['level_01']
                for goid, goobj in self.gosubdag.go2obj.items():
                    if goobj.level == 1:
                        if goid not in goid2color:
                            goid2color[goid] = color
        return goid2color

    def _init_equiv(self):
        """Add equivalent GO IDs to go2color, if necessary."""
        gocolored_all = set(self.go2color)
        go2obj_usr = self.gosubdag.go2obj
        go2color_add = {}
        for gocolored_cur, color in self.go2color.items():
            # Ignore GOs in go2color that are not in the user set
            if gocolored_cur in go2obj_usr:
                goobj = go2obj_usr[gocolored_cur]
                goids_equiv = goobj.alt_ids.union([goobj.id])
                # mrk_alt = "*" if gocolored_cur != goobj.id else ""
                # print("COLORED({}) KEY({}){:1} ALL({})".format(
                #     gocolored_cur, goobj.id, mrk_alt, goids_equiv))
                # Loop through GO IDs which are not colored, but are equivalent to colored GO IDs.
                for goid_add in goids_equiv.difference(gocolored_all):
                    if goid_add in go2color_add:
                        print('**TBD: TWO DIFFERENT COLORS FOR EQUIV GO ID') # pylint: disable=superfluous-parens
                    go2color_add[goid_add] = color
        # print("ADDING {N} GO IDs TO go2color".format(N=len(go2color_add)))
        for goid, color in go2color_add.items():
            self.go2color[goid] = color


# Copyright (C) 2016-2018, DV Klopfenstein, H Tang, All rights reserved.
