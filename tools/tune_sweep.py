#!/usr/bin/env python3
"""
Sweep spotTolerance against ground truth at 20 averaged frames. This is what set the default of
5, replacing 10.

Two-sided objective. Too high and real spots are merged away; too low and spurious noise maxima
both (a) trip the spot proximity filter, culling their real neighbours, and (b) mask out
background pixels, since backgroundMask is built from *all* maxima rather than the surviving
ones. Both ends show up in the background error as well as in the recall.
"""
import json, os, shutil, subprocess, sys
import numpy as np
import tifffile

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from rig import CLASSES, JAVA, classpath          # noqa: E402

FRAMES = 20
MEDIAN_N = 400.0
SEED = 11
THRESHOLD = 6.0       # spotThreshold, so "should be found" means true SNR >= this

SIGMAS = [1.0, 1.5, 2.0, 2.5, 3.0]
DENSITIES = [400, 150]
TOLERANCES = [1, 2, 3, 4, 5, 6, 8, 10, 12, 16, 20]

GEN = os.path.join(os.path.dirname(os.path.abspath(__file__)), "gen_calib.py")
CP = classpath()


# MaximumFinder returns the integer coordinate of the brightest pixel, so localization error
# grows with the PSF width - median 0.82 px at sigma 1 against 1.26 px at sigma 3, worst case
# 3.08. A fixed 1.5 px radius counted 20% of real detections as false positives at sigma 2.5.
# Capped at 3 px, half the closest two truth spots are ever placed, so a detection can still
# only match one of them.
def match_radius(sigma):
    return min(3.0, max(1.5, 1.2 * sigma))


rows = []
for nspots in DENSITIES:
    for sigma in SIGMAS:
        MATCH = match_radius(sigma)
        d = f"tol_{nspots}_{sigma}"
        shutil.rmtree(d, ignore_errors=True)
        subprocess.run(["python3", GEN, str(sigma), str(MEDIAN_N), d, str(SEED)],
                       check=True, stdout=subprocess.DEVNULL,
                       env=dict(os.environ, CALIB_SPOTS=str(nspots), CALIB_FRAMES=str(FRAMES)))

        truth = np.genfromtxt(os.path.join(d, "truth.csv"), delimiter=",", names=True)
        tx, ty, tsnr = truth["cx"], truth["cy"], truth["snr"]
        findable = tsnr >= THRESHOLD
        truth_bg = tifffile.imread(os.path.join(d, "truth_bg.tif")).astype(np.float64)
        sy = np.clip(ty.astype(int), 0, truth_bg.shape[0] - 1)
        sx = np.clip(tx.astype(int), 0, truth_bg.shape[1] - 1)

        for tol in TOLERANCES:
            subprocess.run(
                [JAVA, "-Djava.awt.headless=true", "-Xmx8g",
                 "-cp", os.pathsep.join([".", CLASSES, CP]),
                 "RunTune", d, str(sigma), str(tol), "-1000.0", str(FRAMES), str(THRESHOLD)],
                check=True, capture_output=True, text=True)

            # Kept so the sweep can be re-scored without re-running any Java.
            shutil.copy(os.path.join(d, "sim_spotf_spots.csv"),
                        os.path.join(d, f"spots_tol{tol}.csv"))

            spots = np.genfromtxt(os.path.join(d, "sim_spotf_spots.csv"), delimiter=",",
                                  names=True, dtype=None)
            px = np.atleast_1d(spots["x"]).astype(float)
            py = np.atleast_1d(spots["y"]).astype(float)

            if len(px):
                dist = np.hypot(px[:, None] - tx[None, :], py[:, None] - ty[None, :])
                nearest = dist.argmin(axis=1)
                good = dist[np.arange(len(px)), nearest] <= MATCH
                found = np.zeros(len(tx), dtype=bool)
                found[nearest[good]] = True
                nfp = int((~good).sum())
            else:
                found = np.zeros(len(tx), dtype=bool)
                nfp = 0

            recall = float(found[findable].mean()) if findable.any() else 0.0

            # Background estimate quality where it actually gets subtracted.
            est = tifffile.imread(os.path.join(d, "sim_spotf_bg_smooth.tif")).astype(np.float64)
            err = est[sy, sx] - truth_bg[sy, sx]

            # How much of the field the mask left usable.
            masks = tifffile.imread(os.path.join(d, "sim_spotf_masks.tif")).astype(np.float64)
            trusted = float(((masks[0] > 0) & (masks[1] > 0)).mean())

            rows.append({"spots": nspots, "sigma": sigma, "tol": tol,
                         "final": int(len(px)), "truthFindable": int(findable.sum()),
                         "recall": recall, "fp": nfp, "trusted": trusted,
                         "bgBias": float(err.mean()), "bgRms": float(np.sqrt((err ** 2).mean()))})
            r = rows[-1]
            print(f"  {nspots:4d}sp s{sigma:3.1f} tol{tol:3d}  final {r['final']:4d}"
                  f"  recall {recall:5.3f}  fp {nfp:4d}  trusted {trusted:5.3f}"
                  f"  bgBias {r['bgBias']:+7.3f}  bgRms {r['bgRms']:6.3f}", flush=True)

json.dump(rows, open("tune_sweep.json", "w"), indent=1)
print("\nwrote tune_sweep.json - run analyze_tune.py to summarize")
