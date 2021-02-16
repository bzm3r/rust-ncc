def default_ani_opts():
    return AniOpts(50, False, False, False, True)


class AniOpts:
    def __init__(self, rgtp_scale, label_verts, label_cells,
                 show_trails, follow_group):
        self.rgtp_scale = rgtp_scale
        self.label_verts = label_verts
        self.label_cells = label_cells
        self.show_trails = show_trails
        self.follow_group = follow_group

    def description(self):
        ds = "_RS={}".format(self.rgtp_scale) + \
             "_SV={}".format(self.label_verts) + \
             "_SC={}".format(self.label_cells) + \
             "_ST={}".format(self.show_trails) + \
             "_FG={}".format(self.follow_group)
        return ds
