#!/usr/bin/env python3
"""
Compare shape based doublet rejectors against the prominence filter. This is the evidence behind
docs/adr/0001-no-shape-based-doublet-rejection.md - read that before drawing conclusions here.

All computed offline from foreground = analysis image - background estimate, on the same
labelled fields, scored the same way: best achievable FRR(singles) + FAR(bad), where 1.000 is
what keeping every spot scores.

  ellipticity  |e| from Gaussian weighted second moments. A pair separated by d along one axis
               has moment ratio 1 + d^2/4sigma^2 exactly. Note |e| RAW loses to prominence and
               only wins once divided by its 1/SNR noise floor, because |e| for a single is
               positive definite noise.
  size         the moment trace, which goes as 1 + d^2/8sigma^2
  chi2         reduced residual of a single Gaussian fit. Looks like much the strongest
               aggregate detector and is substantially a brightness cut - check the reported
               correlation against SNR among true singles before believing it.

Run gen_shape.py at constant SNR first, and again with SHAPE_SINGLES/DOUBLETS/AGGREGATES set for
a dirty field: the clean and dirty answers differ, and that difference is the ADR's conclusion.
"""
import math
import numpy as np
import tifffile

SIGMAS = [1.0, 1.5, 2.0, 2.5, 3.0]
MATCH = 2.0
BLACK, GAIN, FRAMES = 5, 1.0, 20


def window(sigma):
    reach = int(math.ceil(2.5 * sigma))
    dx, dy = np.meshgrid(np.arange(-reach, reach + 1), np.arange(-reach, reach + 1))
    return dx.ravel(), dy.ravel()


def patches(img, xs, ys, dx, dy):
    h, w = img.shape
    px = xs[:, None] + dx[None, :]
    py = ys[:, None] + dy[None, :]
    inside = (px >= 0) & (px < w) & (py >= 0) & (py < h)
    return np.where(inside, img[np.clip(py, 0, h - 1), np.clip(px, 0, w - 1)], 0.0)


def measure(img, xs, ys, sigma, noise):
    """Weighted moments, then a linear single Gaussian fit at the measured centroid."""
    dx, dy = window(sigma)
    I = patches(img, xs, ys, dx, dy)

    w = np.exp(-(dx ** 2 + dy ** 2) / (2.0 * sigma * sigma))[None, :]
    wi = w * I
    norm = wi.sum(axis=1)
    norm = np.where(np.abs(norm) < 1e-9, 1e-9, norm)

    cx = (wi * dx[None, :]).sum(axis=1) / norm
    cy = (wi * dy[None, :]).sum(axis=1) / norm
    ux = dx[None, :] - cx[:, None]
    uy = dy[None, :] - cy[:, None]
    mxx = (wi * ux * ux).sum(axis=1) / norm
    myy = (wi * uy * uy).sum(axis=1) / norm
    mxy = (wi * ux * uy).sum(axis=1) / norm

    trace = mxx + myy
    safe = np.where(np.abs(trace) < 1e-9, 1e-9, trace)
    ellip = np.sqrt((mxx - myy) ** 2 + (2.0 * mxy) ** 2) / np.abs(safe)

    g = np.exp(-((ux ** 2 + uy ** 2) / (2.0 * sigma * sigma)))
    one = np.ones_like(g)
    sgg = (g * g).sum(axis=1)
    sg1 = g.sum(axis=1)
    s11 = one.sum(axis=1)
    sgi = (g * I).sum(axis=1)
    s1i = I.sum(axis=1)
    det = sgg * s11 - sg1 * sg1
    det = np.where(np.abs(det) < 1e-9, 1e-9, det)
    amp = (sgi * s11 - s1i * sg1) / det
    off = (sgg * s1i - sg1 * sgi) / det

    resid = I - (amp[:, None] * g + off[:, None])
    var = np.maximum(noise[:, None] ** 2, 1e-6)
    chi2 = (resid ** 2 / var).sum(axis=1) / max(1, I.shape[1] - 2)
    return ellip, trace, chi2


def load(tag, sigma):
    d = f"{tag}_{sigma}"
    bg = tifffile.imread(f"{d}/sim_spotf_bg_smooth.tif").astype(np.float64)
    fg = tifffile.imread(f"{d}/sim_spotf_qc_image.tif").astype(np.float64) - bg
    sp = np.genfromtxt(f"{d}/sim_spotf_spots.csv", delimiter=",", names=True)
    xs = np.atleast_1d(sp["x"]).astype(int)
    ys = np.atleast_1d(sp["y"]).astype(int)
    snr = np.atleast_1d(sp["snr"]).astype(float)
    prom = np.atleast_1d(sp["prominence"]).astype(float)

    tr = np.genfromtxt(f"{d}/truth.csv", delimiter=",", names=True, dtype=None, encoding="utf-8")
    dist = np.hypot(xs[:, None] - tr["cx"][None, :], ys[:, None] - tr["cy"][None, :])
    near = dist.argmin(axis=1)
    ok = dist[np.arange(len(xs)), near] <= MATCH
    kind = np.array([str(k) for k in tr["kind"]])[near]
    sep = (tr["sep"] / sigma)[near]

    noiseImg = np.sqrt(np.maximum(0.0, bg - 2 * BLACK) / (GAIN * FRAMES))
    noise = noiseImg[np.clip(ys, 0, bg.shape[0] - 1), np.clip(xs, 0, bg.shape[1] - 1)]
    ellip, trace, chi2 = measure(fg, xs, ys, sigma, noise)
    return {"kind": kind[ok], "sep": sep[ok], "snr": snr[ok],
            "prominence": prom[ok], "ellip": ellip[ok],
            "ellipN": (ellip * snr)[ok], "size": (trace / np.median(trace))[ok],
            "chi2": chi2[ok]}


def best(values, single, bad, low_is_bad):
    b = 1.0
    for t in np.unique(np.percentile(values, np.arange(0, 100.01, 0.5))):
        if low_is_bad:
            frr = (single & (values < t)).sum() / max(1, single.sum())
            far = (bad & (values >= t)).sum() / max(1, bad.sum())
        else:
            frr = (single & (values > t)).sum() / max(1, single.sum())
            far = (bad & (values <= t)).sum() / max(1, bad.sum())
        b = min(b, frr + far)
    return b


STATS = [("prominence", True), ("ellip", False), ("ellipN", False),
         ("size", False), ("chi2", False)]

print("Best achievable FRR+FAR, clean fields (1.000 = no better than keeping everything)\n")
print("  statistic     " + "".join(f"{f's{s}':>9}" for s in SIGMAS) + "     mean")
for name, low_bad in STATS:
    row = []
    for sigma in SIGMAS:
        d = load("shape", sigma)
        row.append(best(d[name], d["kind"] == "single", d["kind"] != "single", low_bad))
    print(f"  {name:<14}" + "".join(f"{v:>9.3f}" for v in row) + f"   {np.mean(row):>6.3f}")

print("\nAgainst doublets only:")
for name, low_bad in STATS:
    row = []
    for sigma in SIGMAS:
        d = load("shape", sigma)
        single, dbl = d["kind"] == "single", d["kind"] == "doublet"
        keep = single | dbl
        row.append(best(d[name][keep], single[keep], dbl[keep], low_bad))
    print(f"  {name:<14}" + "".join(f"{v:>9.3f}" for v in row) + f"   {np.mean(row):>6.3f}")

print("\nIs chi2 a brightness cut? r(log snr, log chi2) among TRUE SINGLES:")
for sigma in SIGMAS:
    d = load("shape", sigma)
    m = d["kind"] == "single"
    r = np.corrcoef(np.log(np.maximum(d["snr"][m], 1e-3)),
                    np.log(np.maximum(d["chi2"][m], 1e-9)))[0, 1]
    print(f"  sigma {sigma}: {r:+.3f}")

print("\n|e| x SNR response by separation, as a multiple of the single median:")
for sigma in (1.0, 2.0, 3.0):
    d = load("shape", sigma)
    base = np.median(d["ellipN"][d["kind"] == "single"])
    line = []
    for s in sorted(set(d["sep"][d["kind"] == "doublet"])):
        sel = (d["kind"] == "doublet") & (d["sep"] == s)
        if sel.sum() >= 2:
            line.append(f"{s:.2f}s:{np.median(d['ellipN'][sel]) / base:.2f}x")
    print(f"  sigma {sigma}: " + "  ".join(line))
