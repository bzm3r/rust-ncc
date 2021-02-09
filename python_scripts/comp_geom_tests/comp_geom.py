import copy

import numpy as np
import numpy.random as npr
import matplotlib.pyplot as plt

seed = npr.choice(np.arange(0, 10000))
print("seed: ", seed)
npr.seed(seed)


def is_left(a, b, p):
    ab = b - a
    bp = p - b
    area = ab[0] * bp[1] - ab[1] * bp[0]
    if area > 0:
        return 1
    elif area < 0:
        return -1
    else:
        return 0


def is_point_in_poly(p, poly):
    wn = 0
    for vi in range(len(poly)):
        q = poly[vi]
        r = poly[(vi + 1) % len(poly)]
        if (q[1] < p[1] and r[1] < p[1]) or (q[1] > p[1] and r[1] > p[1]) or (
                q[0] < p[0] and r[0] < p[0]):
            continue
        else:
            wn += is_left(q, r, p)

    return wn > 0


def random_to(mag, a, b, left=1):
    ab = b - a
    ab_perp = left * -1 * np.array([ab[1], -1 * ab[0]])
    mag = (0.1 * npr.rand() + 0.95) * mag
    r = 0.5 * (a + b) + mag * unitize(
        (npr.rand() + 0.1) * ab + (npr.rand() + 0.1) * ab_perp)
    return r


def unitize(a):
    return a / np.linalg.norm(a)


def random_poly(edge_mag, nverts):
    # initialize triangle
    polys = []
    poly = [npr.rand(2) for i in range(2)]
    poly.append(random_to(edge_mag, poly[0], poly[1]))
    polys.append(copy.deepcopy(np.array(poly)))

    while len(poly) < nverts:
        okay_vis = []
        if len(poly) == 3:
            okay_vis = [n for n in range(0, 3)]
        else:
            for vi, v in enumerate(poly):
                u = poly[(vi - 1) % len(poly)]
                w = poly[(vi + 1) % len(poly)]
                x = 0.5 * ((v - u) + (w - v)) + u
                if is_point_in_poly(x, poly):
                    okay_vis.append(vi)
        if len(okay_vis) < 1:
            return polys
        vi = npr.choice(okay_vis)
        a, b = poly[vi], poly[(vi + 1) % len(poly)]
        rtr = random_to(edge_mag, a, b, left=-1)
        poly.insert((vi + 1) % len(poly), rtr)
        vi += 1
        polys.append(copy.deepcopy(np.array(poly)))

    return polys


def star_poly(nverts, variation):
    d_theta = 2 * np.pi / nverts
    rs = (1.0 - variation) * npr.rand(nverts) + 2 * variation
    rs = nverts * rs / np.sum(rs)
    return np.array(
        [[r * np.cos(k * d_theta), r * np.sin(k * d_theta)] for k, r in
         enumerate(rs)])


def plot_poly(mpl_ax, poly):
    mpl_ax.plot(list(poly[:, 0]) + [poly[0, 0]],
                list(poly[:, 1]) + [poly[0, 1]])


def calc_bbox(cell):
    return (np.min(cell[:, 0]), np.max(cell[:, 0]), np.min(cell[:, 1]),
            np.max(cell[:, 1]))


def do_bboxes_intersect(bb0, bb1):
    ax_min, ax_max, ay_min, ay_max = bb0
    bx_min, bx_max, by_min, by_max = bb1

    if ax_min > bx_max:
        # print("ax_min > bx_max")
        return False
    elif ax_max < bx_min:
        # print("ax_max < bx_min")
        return False
    elif ay_max < by_min:
        # print("ay_max < by_min")
        return False
    elif ay_min > by_max:
        # print("ay_min < by_max")
        return False
    else:
        return True


def gen_random_star_polys(N, nverts, variation):
    polys = []
    bboxes = []
    while len(polys) < N:
        sp = star_poly(nverts, variation) + 10 * npr.rand(2)
        sp_bb = calc_bbox(sp)
        if not np.any([do_bboxes_intersect(sp_bb, bb) for bb in bboxes]):
            polys.append(sp)
            bboxes.append(sp_bb)
    return np.array(polys)


def calc_centroid(poly):
    return np.array([np.average(poly[:, 0]), np.average(poly[:, 1])])


def skip(alpha, a, b, iota, i, j):
    if alpha == a and (iota == i or iota == (i - 1) % nverts):
        return True
    if alpha == b and (iota == j or iota == (j - 1) % nverts):
        return True
    return False


def lsegs_intersect(a, b, c, d):
    if is_left(a, b, c) != is_left(a, b, d) and is_left(c, d, a) != is_left(c,
                                                                            d,
                                                                            b):
        return True
    else:
        return False


def lseg_self_intersects_poly(a, b, poly, source_ix):
    nverts = len(poly)
    for ix in range(nverts):
        if ix != source_ix and ix != (source_ix - 1) % nverts:
            alpha = poly[ix]
            beta = poly[(ix + 1) % nverts]
            if lsegs_intersect(a, b, alpha, beta):
                return True

    return False


def lseg_intersects_poly(a, b, poly):
    nverts = len(poly)
    for ix in range(nverts):
        alpha = poly[ix]
        beta = poly[(ix + 1) % nverts]
        if lsegs_intersect(a, b, alpha, beta):
            return True

    return False

# def paint():
#     global fig
#     global ax
#     global ix
#     global ns
#     ax.cla()
#     ax.set_aspect('equal')
#     ax.set_xlim(min_n, max_n)
#     ax.set_ylim(min_n, max_n)
#     ax.set_title("poly: {}".format(ix))
#     plot_poly(ax, polys[ix])
#
#
# def on_press(event):
#     global fig, ix, num_polys
#     if event.key == "x":
#         ix = (ix + 1) % num_polys
#     if event.key == "z":
#         ix = (ix - 1) % num_polys
#     paint()
#     fig.canvas.draw()
#
#
# polys = random_poly(1.0, 8)
# num_polys = len(polys)
# ix = 0
# fig, ax = plt.subplots()
# ns = []
# for poly in polys:
#     ns += list(poly.flatten())
# max_n = np.max(ns)
# min_n = np.min(ns)
# fig.canvas.mpl_connect('key_press_event', on_press)
#
fig, ax = plt.subplots()
ax.set_aspect("equal")
nverts = 16
num_polys = 10
polys = gen_random_star_polys(num_polys, nverts, 0.0)
n_min = np.min(polys.flatten()) * 1.1
n_max = np.max(polys.flatten()) * 1.1
for poly in polys:
    plot_poly(ax, poly)
lsegs = []
ok_lsegs = []
for n in range(num_polys):
    for i in range(nverts):
        for m in range(n + 1, num_polys):
            for j in range(nverts):
                lsegs.append((n, i, m, j))
for lseg in lsegs:
    n, i, m, j = lseg
    poly_a = polys[n]
    poly_b = polys[m]
    a, b = poly_a[i], poly_b[j]

    los_penalty = 0.0
    for pi, poly in enumerate(polys):
        if pi == n or pi == m:
            for ix in range(0, nverts):
                if (lseg_self_intersects_poly(a, b, poly_a, i) or lseg_self_intersects_poly(a, b, poly_b, j)):
                    los_penalty += 10000
                    break
        else:
            for ix in range(0, nverts):
                if lseg_intersects_poly(a, b, poly):
                    los_penalty += 1
                    break

        ax.plot([a[0], b[0]], [a[1], b[1]], color="k", linewidth=0.25 * (1 / (1 + los_penalty)))
