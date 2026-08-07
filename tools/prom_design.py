#!/usr/bin/env python3
"""
Settle the prominence definition and its threshold. This is what produced the current filter.

Every candidate is normalized by its own noiseless ideal - the same ring statistic applied to a
unit height Gaussian of that sigma over the same pixel offsets - so a perfect isolated spot
scores 1.0 at every spot size and the threshold reads as a fraction of ideal.

Candidates:
  thick/max   the shipped ring before this, (srad-1)^2 <= r^2 <= (srad+1)^2 summarised by its
              brightest pixel
  thin /max   a thin annulus at 2 sigma, still the brightest pixel
  thin /p90   the annulus summarised by its 90th percentile   <- what won
  *  smoothed  the same on the smoothed foreground rather than the raw one

Run gen_shape.py at constant SNR first - `N` proportional to sigma - or the large sigma fields
sit at the detection limit and the comparison is meaningless.
"""
import math
import numpy as np
import tifffile

SIGMAS = [1.0, 1.5, 2.0, 2.5, 3.0]
MATCH = 2.0
HALF = 0.6          # annulus half width, pixels
BLACK, GAIN, FRAMES = 5, 1.0, 20
Z = 1.2816          # 90th percentile of a standard normal


def annulus(sigma):
    out = []
    reach = int(math.ceil(2.0 * sigma + HALF)) + 1
    for rx in range(-reach, reach + 1):
        for ry in range(-reach, reach + 1):
            r = math.hypot(rx, ry)
            if abs(r - 2.0 * sigma) <= HALF:
                out.append((rx, ry, r))
    return out


def thick_ring(sigma):
    srad = int(round(2.0 * sigma))
    return [(rx, ry, math.hypot(rx, ry))
            for rx in range(-srad, srad + 1) for ry in range(-srad, srad + 1)
            if (srad - 1) ** 2 <= rx * rx + ry * ry <= (srad + 1) ** 2]


def ideal_ratio(offsets, sigma, how):
    """Ring statistic for a noiseless unit height Gaussian, inverted."""
    profile = np.array([math.exp(-(r * r) / (2.0 * sigma * sigma)) for _, _, r in offsets])
    ref = profile.max() if how == "max" else np.percentile(profile, 90)
    return 1.0 / ref


def value_at(img, x, y):
    h, w = img.shape
    inside = (x >= 0) & (x < w) & (y >= 0) & (y < h)
    return np.where(inside, img[np.clip(y, 0, h - 1), np.clip(x, 0, w - 1)], 0.0)


def prominence(img, xs, ys, offsets, how, ideal, noise=None):
    height = value_at(img, xs, ys)
    ring = np.stack([value_at(img, xs + rx, ys + ry) for rx, ry, _ in offsets], axis=1)
    ref = ring.max(axis=1) if how == "max" else np.percentile(ring, 90, axis=1)

    # The ring pixels are background subtracted, so their noise is zero mean and the 90th
    # percentile of it sits Z sigmas above the true ring level. Left in, that reads the ring high
    # and the prominence low, and by more at larger spot sizes - which is what stops the
    # normalization above from actually equalizing the scores.
    if noise is not None:
        ref = ref - Z * noise

    # A ring at or below zero is a ring with nothing in it, so there is no neighbour.
    raw = np.where(ref > 0.0, height / np.maximum(ref, 1.0e-9), np.inf)
    return np.minimum(raw / ideal, 10.0)


def load(tag, sigma):
    d = f"{tag}_{sigma}"
    bg = tifffile.imread(f"{d}/sim_spotf_bg_smooth.tif").astype(np.float64)
    raw = tifffile.imread(f"{d}/sim_spotf_qc_image.tif").astype(np.float64) - bg
    smooth = tifffile.imread(f"{d}/sim_spotf_fg_smooth.tif").astype(np.float64)
    noiseImg = np.sqrt(np.maximum(0.0, bg - 2 * BLACK) / (GAIN * FRAMES))

    sp = np.genfromtxt(f"{d}/sim_spotf_spots.csv", delimiter=",", names=True)
    xs = np.atleast_1d(sp["x"]).astype(int)
    ys = np.atleast_1d(sp["y"]).astype(int)

    tr = np.genfromtxt(f"{d}/truth.csv", delimiter=",", names=True, dtype=None, encoding="utf-8")
    dist = np.hypot(xs[:, None] - tr["cx"][None, :], ys[:, None] - tr["cy"][None, :])
    near = dist.argmin(axis=1)
    ok = dist[np.arange(len(xs)), near] <= MATCH
    kind = np.array([str(k) for k in tr["kind"]])[near]
    noise = noiseImg[np.clip(ys, 0, bg.shape[0] - 1), np.clip(xs, 0, bg.shape[1] - 1)]
    return raw, smooth, xs[ok], ys[ok], kind[ok], noise[ok]


def best_rate_sum(values, single, bad):
    best = 1.0
    for t in np.unique(np.percentile(values, np.arange(0, 100.01, 0.5))):
        frr = (single & (values < t)).sum() / max(1, single.sum())
        far = (bad & (values >= t)).sum() / max(1, bad.sum())
        best = min(best, frr + far)
    return best


CANDIDATES = [
    ("thick/max  raw   (old)", "thick", "max", "raw"),
    ("thin /max  raw", "thin", "max", "raw"),
    ("thin /p90  raw", "thin", "p90", "raw"),
    ("thin /max  smoothed", "thin", "max", "smooth"),
    ("thin /p90  smoothed", "thin", "p90", "smooth"),
]

print("Best achievable FRR+FAR per definition (1.000 = no better than keeping every spot)\n")
print(f"  {'definition':<26}" + "".join(f"{s:>8}" for s in SIGMAS) + "     mean")
for label, ring, how, src in CANDIDATES:
    row = []
    for sigma in SIGMAS:
        raw, smooth, xs, ys, kind, noise = load("shape", sigma)
        offs = annulus(sigma) if ring == "thin" else thick_ring(sigma)
        img = raw if src == "raw" else smooth
        v = prominence(img, xs, ys, offs, how, ideal_ratio(offs, sigma, how))
        row.append(best_rate_sum(v, kind == "single", kind != "single"))
    print(f"  {label:<26}" + "".join(f"{x:>8.3f}" for x in row) + f"   {np.mean(row):>6.3f}")

print("\nShipped definition - thin annulus, 90th percentile, noise bias removed.")
print("Median score by object, and the error count at each threshold.\n")
data = {}
for sigma in SIGMAS:
    raw, smooth, xs, ys, kind, noise = load("shape", sigma)
    offs = annulus(sigma)
    v = prominence(raw, xs, ys, offs, "p90", ideal_ratio(offs, sigma, "p90"), noise)
    data[sigma] = (v, kind)
    print(f"  sigma {sigma}: single {np.median(v[kind == 'single']):.3f}"
          f"   doublet {np.median(v[kind == 'doublet']):.3f}"
          f"   aggregate {np.median(v[kind == 'aggregate']):.3f}")

print("\n  thr   " + "".join(f"{f's{s}':>8}" for s in SIGMAS) + "    total   singles kept")
for t in [0.0, 0.3, 0.4, 0.5, 0.6, 0.65, 0.7, 0.8, 0.9]:
    cells, total, kept, allsingle = [], 0, 0, 0
    for sigma in SIGMAS:
        v, kind = data[sigma]
        single, bad = kind == "single", kind != "single"
        lost = int((single & (v < t)).sum())
        err = lost + int((bad & (v >= t)).sum())
        total += err
        kept += int(single.sum()) - lost
        allsingle += int(single.sum())
        cells.append(f"{err:>8d}")
    print(f"  {t:4.2f}  " + "".join(cells) + f"{total:>9d}   {100.0 * kept / allsingle:5.1f}%")
print("\nThreshold 0 is the filter off. The default of 0.4 is the largest value that costs no "
      "true singles at all up to sigma 2.")
