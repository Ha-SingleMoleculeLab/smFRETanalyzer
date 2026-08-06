# ADR-0001: No shape based doublet rejection in the spot finder

- **Status:** accepted
- **Date:** 2026-08-06
- **Affects:** `smFRETSpotFinder.spotFilterProminence`

## Context

Re-deriving `spotProminence` established that the prominence filter cannot see the thing it
exists to catch. Two molecules closer than about one sigma are simply a brighter point source,
so a tight pair scores *higher* than a single, not lower. That left an open question: is there a
better measurement for pairs separated by roughly 1 to 3 sigma - close enough that
`MaximumFinder` returns one maximum, so the spot proximity filter never sees them, but far
enough that the combined image is not a single Gaussian?

The physics is favourable on paper. Two equal molecules separated by `d` along one axis give a
second moment of `sigma^2 + d^2/4` along that axis against `sigma^2` across it, so the moment
ratio is exactly `1 + d^2/4 sigma^2` - 2.00 at a separation of 2 sigma, an axis ratio of 1.41.
That is a far larger signal than anything a brightness ring carries. Unequal pairs are weaker,
`sigma^2 + f(1-f) d^2`, so a 20:80 pair gives 0.16 d^2 rather than 0.25 d^2.

## What was measured

Four statistics, all computed offline from `foreground = analysis image - background estimate`
on the labelled simulated fields (singles, doublets at 0.25 to 3 sigma, aggregates), scored the
same way: best achievable `FRR(singles) + FAR(bad)` over a threshold scan, where 1.000 is what
keeping every spot scores.

- **prominence** - what the plugin ships, a thin ring at 2 sigma summarised by its 90th
  percentile, noise bias removed, over the noiseless ideal.
- **ellipticity** - `|e|` from Gaussian weighted second moments, raw and multiplied by SNR.
- **size** - the moment trace `Mxx + Myy`, which goes as `1 + d^2/8 sigma^2`.
- **chi2** - reduced residual of a single Gaussian fit. Sigma is known and the centroid comes
  from the moments, so amplitude and a local offset are a 2x2 linear solve, not a nonlinear fit.

### Against doublets, clean fields

| statistic | s1.0 | s1.5 | s2.0 | s2.5 | s3.0 | mean |
|---|---|---|---|---|---|---|
| prominence (shipped) | 0.396 | 0.457 | 0.576 | **0.498** | 0.765 | 0.538 |
| **\|e\| x SNR** | **0.279** | **0.317** | **0.427** | 0.646 | **0.698** | **0.473** |
| \|e\| raw | 0.542 | 0.494 | 0.602 | 0.872 | 0.894 | 0.681 |
| chi2 | 0.311 | 0.355 | 0.679 | 0.715 | 0.759 | 0.564 |

`|e| x SNR` is the best of the four, winning at four spot sizes of five.

### Two results that are easy to get wrong

**Raw ellipticity loses to prominence, and only wins after being divided by its noise floor.**
`|e|` for a single is positive definite noise whose floor scales as 1/SNR, so pooling spots of
different brightness smears the distribution and destroys the separation the median shift
suggests. Normalizing moves it from 0.681 to 0.473. Anyone reaching for moments here must
normalize; the mature version of this is astronomy's adaptive moment machinery.

**chi2 is substantially a brightness cut.** Among *true singles*, `r(log SNR, log chi2)` is
**+0.67** at sigma 1, +0.37 at sigma 2, +0.44 at sigma 3. Its striking aggregate response - 55x
the single median at sigma 1 - is largely "aggregates are bright" rather than "aggregates are
the wrong shape". This is the same defect the old prominence had with the sign reversed, and it
is why chi2 is not recommended despite scoring well on paper.

### Where the signal actually is

Median `|e| x SNR` as a multiple of the single median, at sigma 2.0:

| separation | 0.25s | 0.50s | 0.75s | 1.00s | 1.25s | 1.50s | 2.00s | aggregate |
|---|---|---|---|---|---|---|---|---|
| response | 0.93x | 0.86x | 2.11x | 3.35x | 3.17x | 3.61x | 4.65x | 13.71x |

The usable band starts at **0.75 sigma**, not 1 sigma. Below 0.5 sigma there is no signal at
all - the response is at or below the single median, which is what the physics predicts and no
shape measure can escape.

## Decision

**Do not add a shape based filter.** Keep the prominence filter as re-derived, and treat the
sub-1-sigma band as out of reach of the spot finder.

The reason is not that the statistic is bad - it is measurably better than what ships. It is
that a *constant* threshold, which is what a default has to be, buys almost nothing on the
fields this pipeline actually sees:

| field | filter off | best constant | gain |
|---|---|---|---|
| clean (345 : 50 : 5) | 208 errors | `\|e\| x SNR <= 2.0` -> 198 | **5%** |
| dirty (200 : 100 : 10) | 233 errors | `\|e\| x SNR <= 0.75` -> 144 | **38%** |

The two optima are far apart and the wrong one does real damage: 0.75 on a clean field gives
290 errors against 208 for no filter at all. This is structural rather than a property of this
particular statistic - when good spots outnumber bad roughly 7 to 1, nearly any filter loses
more singles than it saves bad objects. The example fields are clean, so the 5% column is the
one that applies.

## Consequences

- Pairs closer than about 0.75 sigma remain undetected by the spot finder, and pairs from 0.75
  to 3 sigma are only weakly filtered. This is accepted, not overlooked.
- **The real answer for that band is two step photobleaching**, which works at any separation
  including the sub-0.75-sigma range no shape measure can reach: a single donor bleaches in one
  step, a pair in two. It needs the whole trace, so it belongs in `smFRETAnalyzer` or the
  histogram viewer rather than the spot finder. An integrated intensity cut is a weaker version
  of the same idea and is already partly available through the histogram's intensity range.
- If a shape filter is ever added, it should be **alongside** prominence rather than replacing
  it. Against aggregates the two score 0.092 and 0.084, essentially tied, but for different
  reasons, and the simulation flatters ellipticity there: the aggregates are random clusters of
  4 to 10 molecules and so are anisotropic by chance. A genuinely symmetric aggregate would be
  invisible to ellipticity and visible to prominence.

## What would reopen this

- Real fields substantially dirtier than the example movies. The 38% figure is real; it just
  does not apply to clean data. If the bad fraction approaches a third, add the filter.
- A per spot size threshold rather than a constant. Optimized per sigma the gains are 17% for
  prominence and 41% for `|e| x SNR` on clean fields; the whole loss is in forcing one number.
- Two step photobleaching being implemented, which would change what the spot finder needs to
  catch at all.

## Caveats on these numbers

The aggregate samples are 4 to 5 objects per field, so the aggregate columns are indicative
only. "Best achievable" is maximized over roughly 200 thresholds on the same data and is
therefore optimistically biased; differences under about 0.05 should not be trusted.

Both measurements depend on two methodology points recorded in `CLAUDE.md`: simulated fields
must hold SNR constant rather than brightness, and the detection to truth match radius must
scale with sigma.
