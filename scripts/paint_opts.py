def default_ani_opts():
    return AniOpts(50, False, False, False, True)


class AniOpts:
    def __init__(self, arrow_scale, label_verts, label_cells,
                 show_trails, show_group_trail, follow_group):
        self.arrow_scale = arrow_scale
        self.label_verts = label_verts
        self.label_cells = label_cells
        self.show_trails = show_trails
        self.follow_group = follow_group
        self.show_group_trail = show_group_trail
        self.poly_line_style = "-"

    def description(self):
        ds = "_RS={}".format(self.arrow_scale) + \
             "_SV={}".format(self.label_verts) + \
             "_SC={}".format(self.label_cells) + \
             "_ST={}".format(self.show_trails) + \
             "_FG={}".format(self.follow_group) + \
             "_SG={}".format(self.show_group_trail)
        return ds
