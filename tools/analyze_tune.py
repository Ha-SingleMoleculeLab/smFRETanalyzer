#!/usr/bin/env python3
"""Summarize tune_sweep.json: the plateau per condition, and whether one constant fits them all."""
import json
import numpy as np

with open("tune_sweep.json") as fh:
    rows = json.load(fh)
sigmas = sorted(set(r["sigma"] for r in rows))
densities = sorted(set(r["spots"] for r in rows), reverse=True)
tols = sorted(set(r["tol"] for r in rows))


def table(title, key, fmt):
    print(f"\n{title}\n")
    print("          " + "".join(f"{t:>7d}" for t in tols))
    for nsp in densities:
        for s in sigmas:
            sub = {r["tol"]: r for r in rows if r["spots"] == nsp and r["sigma"] == s}
            print(f"{nsp:4d} s{s:3.1f} " + "".join(format(sub[t][key], fmt) for t in tols))


table("Recall (fraction of truth spots with true SNR >= 6 that survive every filter)",
      "recall", ">7.3f")
table("False positives (detections with no truth spot inside the match radius)", "fp", ">7d")
table("Background rms error at spot locations, ADU", "bgRms", ">7.3f")
table("Fraction of the field the background estimator was allowed to use", "trusted", ">7.3f")

plateaus = {}
for nsp in densities:
    for s in sigmas:
        sub = {r["tol"]: r for r in rows if r["spots"] == nsp and r["sigma"] == s}
        best = max(r["recall"] for r in sub.values())
        plateaus[(nsp, s)] = [t for t in tols if sub[t]["recall"] >= 0.98 * best]

common = set(tols)
for p in plateaus.values():
    common &= set(p)
print("\nTolerances inside every condition's plateau:", sorted(common) if common else "NONE")

print("\nMean recall as a fraction of the best available in that condition:")
scores = []
for t in tols:
    rel = []
    for nsp in densities:
        for s in sigmas:
            sub = {r["tol"]: r for r in rows if r["spots"] == nsp and r["sigma"] == s}
            best = max(r["recall"] for r in sub.values())
            rel.append(sub[t]["recall"] / best if best > 0 else 0.0)
    scores.append((t, float(np.mean(rel)), float(np.min(rel))))
    print(f"  tol {t:3d}   mean {scores[-1][1]:.3f}   worst condition {scores[-1][2]:.3f}")
print("\nBest by mean:", max(scores, key=lambda x: x[1])[0],
      " best by worst case:", max(scores, key=lambda x: x[2])[0])

print("\nPlateau edges against the 1/sigma^2 peak amplitude argument (400 spot fields):")
for s in sigmas:
    p = plateaus[(400, s)]
    print(f"  sigma {s:3.1f}  plateau {min(p):2d} .. {max(p):2d}   "
          f"peak amplitude of a median spot ~ {400.0 / (2 * np.pi * s * s):5.1f} ADU")
