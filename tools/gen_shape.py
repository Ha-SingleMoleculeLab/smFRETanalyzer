#!/usr/bin/env python3
"""
Simulated movies for calibrating the prominence filter: singles, unresolved doublets and
aggregates on an otherwise clean field.

Same illumination, brightness distribution and noise as gen_calib.py. What is added is a
labelled population of things a spot filter is supposed to reject:

  single     one molecule, should be kept
  doublet    two molecules a controlled separation apart, should be rejected when they merge
             into one maximum - the spot proximity filter cannot see them, it only ever gets
             one detection
  aggregate  several molecules inside about a sigma, so a bright slightly extended blob

The default mix is roughly ten doublets per aggregate on a field that is mostly singles, which
is what the real fields look like. Raise SHAPE_DOUBLETS and SHAPE_AGGREGATES for a dirty field;
the two regimes want very different filter settings and that difference is the main finding in
docs/adr/0001-no-shape-based-doublet-rejection.md.

Usage: gen_shape.py <sigma> <median N> <outdir> <seed>
Environment: SHAPE_SINGLES, SHAPE_DOUBLETS, SHAPE_AGGREGATES, SHAPE_FRAMES.
"""
import json, math, os, sys
import numpy as np
import tifffile

# scipy was pulled in for erf alone. These arrays are a few tens of elements per spot, so a
# vectorized math.erf costs nothing measurable and drops the heaviest dependency. The two erfs
# are not bitwise equal - they differ by up to one ulp, 3.3e-16 - but that is far below the
# Poisson quantization that follows, and swapping one for the other was checked to produce
# byte-identical movies from both generators.
erf = np.vectorize(math.erf, otypes=[np.float64])

W, H = 256, 512
BLACK = 5
BACK = 20
BEAM = 180.0
FLOOR = 0.5
MIN_SEP = 8.0             # between objects, pixels - clean field, objects well separated
SPREAD = 0.5              # sigma of ln(N)

N_SINGLE = int(os.environ.get("SHAPE_SINGLES", "345"))
N_DOUBLET = int(os.environ.get("SHAPE_DOUBLETS", "50"))
N_AGGREGATE = int(os.environ.get("SHAPE_AGGREGATES", "5"))
FRAMES = int(os.environ.get("SHAPE_FRAMES", "20"))

# Separations to spread the doublets over, in units of sigma. Below about 0.5 nothing can
# tell a doublet from a single; above about 3 they resolve and the proximity filter takes over.
SEPARATIONS = [0.25, 0.5, 0.75, 1.0, 1.25, 1.5, 2.0, 2.5, 3.0]


def illumination():
    y, x = np.mgrid[0:H, 0:W].astype(np.float64)
    r2 = (x - W / 2.0) ** 2 + (y - H / 2.0) ** 2
    return FLOOR + (1.0 - FLOOR) * np.exp(-r2 / (2.0 * BEAM * BEAM))


def place(rng, n):
    centres = []
    tries = 0
    while (len(centres) < n) and (tries < n * 600):
        tries += 1
        cx = rng.uniform(14, W - 14)
        cy = rng.uniform(14, H - 14)
        if centres and (np.hypot(*(np.array(centres) - [cx, cy]).T).min() < MIN_SEP):
            continue
        centres.append((cx, cy))
    return np.array(centres)


def render(points, totals, sigma):
    plane = np.zeros((H, W), dtype=np.float64)
    reach = int(np.ceil(5 * sigma)) + 2
    root2 = np.sqrt(2.0) * sigma
    for (cx, cy), total in zip(points, totals):
        x0, x1 = max(0, int(cx) - reach), min(W, int(cx) + reach + 1)
        y0, y1 = max(0, int(cy) - reach), min(H, int(cy) + reach + 1)
        ex = 0.5 * erf((np.arange(x0, x1 + 1) - cx) / root2)
        ey = 0.5 * erf((np.arange(y0, y1 + 1) - cy) / root2)
        plane[y0:y1, x0:x1] += total * np.outer(np.diff(ey), np.diff(ex))
    return plane


def build(sigma, median_n, outdir, seed):
    os.makedirs(outdir, exist_ok=True)
    rng = np.random.default_rng(seed)
    illum = illumination()

    n_objects = N_SINGLE + N_DOUBLET + N_AGGREGATE
    centres = place(rng, n_objects)
    rng.shuffle(centres)

    kinds = (["single"] * N_SINGLE + ["doublet"] * N_DOUBLET + ["aggregate"] * N_AGGREGATE)
    kinds = kinds[:len(centres)]

    points, totals, records = [], [], []
    for (cx, cy), kind in zip(centres, kinds):
        at = illum[int(cy), int(cx)]
        brightness = median_n * np.exp(rng.normal(0.0, SPREAD)) * at

        if kind == "single":
            points.append((cx, cy)); totals.append(brightness)
            records.append((cx, cy, kind, 0.0, brightness, 1))

        elif kind == "doublet":
            sep = SEPARATIONS[len(records) % len(SEPARATIONS)] * sigma
            angle = rng.uniform(0.0, 2.0 * np.pi)
            dx, dy = 0.5 * sep * np.cos(angle), 0.5 * sep * np.sin(angle)
            second = median_n * np.exp(rng.normal(0.0, SPREAD)) * at
            points.append((cx - dx, cy - dy)); totals.append(brightness)
            points.append((cx + dx, cy + dy)); totals.append(second)
            records.append((cx, cy, kind, sep, brightness + second, 2))

        else:
            n_mol = int(rng.integers(4, 11))
            total = 0.0
            for _ in range(n_mol):
                ox, oy = rng.normal(0.0, 0.8 * sigma, 2)
                one = median_n * np.exp(rng.normal(0.0, SPREAD)) * at
                points.append((cx + ox, cy + oy)); totals.append(one)
                total += one
            records.append((cx, cy, kind, 0.0, total, n_mol))

    back = BACK * illum
    rate = render(points, np.array(totals) / 2.0, sigma) + back

    stack = np.empty((FRAMES, 512, 512), dtype=np.uint16)
    for f in range(FRAMES):
        stack[f, :, :W] = rng.poisson(rate) + BLACK
        stack[f, :, W:] = rng.poisson(rate) + BLACK
    tifffile.imwrite(os.path.join(outdir, "sim.tif"), stack)

    pts = [[128.0, 128.0], [64.0, 384.0], [192.0, 384.0], [0.0, 0.0]]
    with open(os.path.join(outdir, "sim_mapping.json"), "w") as fh:
        json.dump({"source points": pts, "target points": pts,
                   "image width": 512, "image height": 512}, fh)
    tifffile.imwrite(os.path.join(outdir, "truth_bg.tif"),
                     (2.0 * (back + BLACK)).astype(np.float32))

    with open(os.path.join(outdir, "truth.csv"), "w") as fh:
        fh.write("cx,cy,kind,sep,total,nmol,snr\n")
        for cx, cy, kind, sep, total, nmol in records:
            b = 2.0 * back[int(cy), int(cx)]
            snr = total / np.sqrt(4.0 * np.pi * sigma * sigma * b)
            fh.write(f"{cx},{cy},{kind},{sep},{total},{nmol},{snr}\n")

    counts = {k: sum(1 for r in records if r[2] == k) for k in ("single", "doublet", "aggregate")}
    with open(os.path.join(outdir, "config.json"), "w") as fh:
        json.dump({"sigma": sigma, "median n": median_n, "frames": FRAMES, **counts}, fh)
    return counts


if __name__ == "__main__":
    c = build(float(sys.argv[1]), float(sys.argv[2]), sys.argv[3], int(sys.argv[4]))
    print(f"{sys.argv[3]}: sigma {sys.argv[1]}, {c['single']} singles, "
          f"{c['doublet']} doublets, {c['aggregate']} aggregates")
