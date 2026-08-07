#!/usr/bin/env python3
"""
Simulated movies for calibrating the background estimator and the spot filters against a known
answer.

Random placement at the density of the real data, log-normal brightness, and a broad Gaussian
beam rather than flat illumination so the sigma 14 background smoothing has something to do.
Spot brightness follows the beam, as real excitation would.

Writes, into <outdir>:
  sim.tif           the movie, both channels side by side, Poisson noise, camera offset added
  sim_mapping.json  an identity mapping, so the split needs no registration
  truth_bg.tif      the true background of the summed image, which is what backgroundEstimate
                    is trying to recover
  truth.csv         cx, cy, total, snr for every spot

Usage: gen_calib.py <sigma> <median N> <outdir> <seed>
Environment: CALIB_SPOTS, CALIB_FRAMES, CALIB_SPREAD.

Note when sweeping over sigma: hold the SNR constant rather than the brightness, or the large
sigma fields sit at the detection limit and every filter looks useless there. See the note in
CLAUDE.md - `N` proportional to sigma keeps the true SNR flat.
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
BLACK = 5                 # camera offset, per channel
BACK = 20                 # peak background above black, per channel, photons/pixel/frame
BEAM = 180.0              # illumination scale, pixels
FLOOR = 0.5               # illumination at the corners, as a fraction of peak
MIN_SEP = 6.0             # closest two spot centres may be, pixels

NSPOTS = int(os.environ.get("CALIB_SPOTS", "400"))
FRAMES = int(os.environ.get("CALIB_FRAMES", "30"))
SPREAD = float(os.environ.get("CALIB_SPREAD", "0.5"))   # sigma of ln(N)


def illumination():
    y, x = np.mgrid[0:H, 0:W].astype(np.float64)
    r2 = (x - W / 2.0) ** 2 + (y - H / 2.0) ** 2
    return FLOOR + (1.0 - FLOOR) * np.exp(-r2 / (2.0 * BEAM * BEAM))


def place(rng):
    centres = []
    tries = 0
    while (len(centres) < NSPOTS) and (tries < NSPOTS * 400):
        tries += 1
        cx = rng.uniform(12, W - 12)
        cy = rng.uniform(12, H - 12)
        if centres:
            if np.hypot(*(np.array(centres) - [cx, cy]).T).min() < MIN_SEP:
                continue
        centres.append((cx, cy))
    return np.array(centres)


def render(centres, totals, sigma):
    plane = np.zeros((H, W), dtype=np.float64)
    reach = int(np.ceil(5 * sigma)) + 2
    root2 = np.sqrt(2.0) * sigma
    for (cx, cy), total in zip(centres, totals):
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
    centres = place(rng)
    at_spot = illum[centres[:, 1].astype(int), centres[:, 0].astype(int)]
    totals = median_n * np.exp(rng.normal(0.0, SPREAD, len(centres))) * at_spot

    back = BACK * illum
    rate = render(centres, totals / 2.0, sigma) + back      # half the total per channel

    stack = np.empty((FRAMES, 512, 512), dtype=np.uint16)
    for f in range(FRAMES):
        stack[f, :, :W] = rng.poisson(rate) + BLACK
        stack[f, :, W:] = rng.poisson(rate) + BLACK
    tifffile.imwrite(os.path.join(outdir, "sim.tif"), stack)

    pts = [[128.0, 128.0], [64.0, 384.0], [192.0, 384.0], [0.0, 0.0]]
    with open(os.path.join(outdir, "sim_mapping.json"), "w") as fh:
        json.dump({"source points": pts, "target points": pts,
                   "image width": 512, "image height": 512}, fh)

    # True background of the summed image, which is what backgroundEstimate estimates.
    tifffile.imwrite(os.path.join(outdir, "truth_bg.tif"),
                     (2.0 * (back + BLACK)).astype(np.float32))

    b_at = 2.0 * back[centres[:, 1].astype(int), centres[:, 0].astype(int)]
    snr = totals / np.sqrt(4.0 * np.pi * sigma * sigma * b_at)
    with open(os.path.join(outdir, "truth.csv"), "w") as fh:
        fh.write("cx,cy,total,snr\n")
        for (cx, cy), t, s in zip(centres, totals, snr):
            fh.write(f"{cx},{cy},{t},{s}\n")

    with open(os.path.join(outdir, "config.json"), "w") as fh:
        json.dump({"sigma": sigma, "median n": median_n, "spots": len(centres),
                   "frames": FRAMES, "snr median": float(np.median(snr))}, fh)
    return len(centres), float(np.median(snr))


if __name__ == "__main__":
    n, s = build(float(sys.argv[1]), float(sys.argv[2]), sys.argv[3], int(sys.argv[4]))
    print(f"{sys.argv[3]}: sigma {sys.argv[1]}, {n} spots, median true SNR {s:.2f}")
