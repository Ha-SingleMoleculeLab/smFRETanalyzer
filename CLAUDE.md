# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project overview

A multipart FIJI/ImageJ plugin (Maven, `pom-scijava` parent) for analyzing two-channel single-molecule FRET (smFRET) microscopy data. It ships as three SciJava `Command` plugins that are meant to be run in sequence, each consuming the previous step's output file(s):

1. **smFRETChannelMapper** (`smFRETChannelMapper.java`) — computes an affine mapping between the donor/target (left half of the image) and acceptor/source (right half) channels using TurboReg landmark registration on an averaged image. Outputs `<name>_mapping.json` (affine transform + expected image size) and a QC overlay TIF.
2. **smFRETSpotFinder** (`smFRETSpotFinder.java`) — loads the mapping JSON, finds candidate spots (local maxima) in the averaged, channel-summed image, then filters them by edge proximity, inter-spot proximity, SNR, and prominence. Outputs `<name>_spotf_finding.json` (parameters + file references), `<name>_spotf_spots.csv` (surviving spot table), `<name>_spotf_qc_image.tif`, and `<name>_spotf_masks.tif` (overlap mask + background mask as a 2-frame stack).
3. **smFRETAnalyzer** (`smFRETAnalyzer.java`) — loads the spot-finder JSON, re-derives the mapping/masks/spot list from it, estimates a per-frame background via temporal boxcar filtering (`Filters3D`), and measures per-spot donor/acceptor time traces. Outputs `<name>.h5` (metadata, spot table, target/source trace matrices) and `<name>.traces` (Taekjip Ha lab binary format).

Each plugin's `run()` method derives a `saveRootName` from its primary input file path (strips the extension) and writes all outputs alongside the input using that root. See `README.md` for the full parameter/output reference exposed to end users in the ImageJ UI.

## Architecture notes

- **Composition, not shared base class.** `smFRETSpotFinder` owns a private `smFRETChannelMapper` instance (`smfcm`) and `smFRETAnalyzer` owns a private `smFRETSpotFinder` instance (`smfsf`). Downstream plugins call into the upstream plugin's public methods (`splitImagePlus`, `loadMappingJSON`, `backgroundEstimate`, etc.) directly on these instances rather than through interfaces. When changing a method signature in one class, check both downstream classes for direct callers.
- **Two-channel image convention.** Every FRET stack is a single image whose left half is the donor/target channel and right half is the acceptor/source channel. `splitImagePlus` (in `smFRETChannelMapper`) is the canonical way to split and optionally TurboReg-transform the source half onto the target half's coordinate frame; both other classes call it rather than re-implementing the split.
- **Mask conventions in `smFRETSpotFinder`:** `overlapMask` marks pixels inside the valid overlap region of both channels (used to drop spots near channel edges); `backgroundMask` marks pixels not within `spotMargin` of any spot (used as the fill target for background inpainting). Both are persisted together as a 2-frame TIF (`_spotf_masks.tif`) and reloaded via `loadMasks()`.
- **Background estimation has two flavors:** `backgroundEstimate(image)` does a full mask-based Gaussian inpaint (repeated `blurGaussian` convolution up to 200 iterations, converging when max pixel delta < 1); `backgroundEstimate(image, fillEstimate)` seeds from a previous frame's estimate first — used by `smFRETAnalyzer.backGroundEstimation()` to amortize inpainting cost across a whole stack (temporal boxcar via `Filters3D.filter`, then per-frame refinement).
- **Spot arrays are `double[][]` shaped `[spot][field]`, and there are two incompatible layouts — check which one you're holding.**
  - *In memory during spot finding:* column 0 is a boolean-as-double "is this spot still good" flag written by each filter stage (`spotFilterWithMask`, `spotFilterSNR`, `spotFilterProminence`), so x/y live at columns 1/2. `getMaxima` returns width 3 (`flag, x, y`); `spotFilterSNR` appends SNR (width 4); `spotFilterProminence` appends prominence (width 5). Filters append a diagnostic column rather than overwriting existing ones.
  - *After a save/load round-trip:* the flag is gone. `saveSpotLocations` writes only good spots starting from index `j+1`, and `loadSpotLocations` reads back into index `j`, so the reloaded array is `x, y, snr, prominence` with **x at column 0** (matching `columnHeaders`).
  - `smFRETAnalyzer` only ever sees the reloaded form — hence `spots[j][0]`/`spots[j][1]` as x/y in `measureTimeTraces` — while everything inside `smFRETSpotFinder` uses the flag-prefixed form. The round trip only lines up because exactly two widening filters run and `columnHeaders.size() == 4`; adding a filter stage means updating `columnHeaders` in lockstep.
- **JSON as the inter-plugin contract.** Jackson `ObjectMapper` reads/writes plain `Map<String,Object>` (via `HashMap.class`), not typed DTOs. Downstream plugins read specific keys by string (e.g. `smFRETAnalyzer.run()` pulls `"root name"`, `"image name"`, `"mapping file"`, `"masks file"`, `"spots file"`, `"spot sigma"`, `"camera black"`, `"camera gain"` out of the spot-finder JSON) and must stay in sync with the keys the upstream plugin writes.
- **TurboReg integration** happens via reflection (`IJ.runPlugIn("TurboReg_", options)` + `getClass().getMethod(...)`) since TurboReg isn't a compile-time Maven dependency — it must be installed separately into the Fiji/ImageJ plugins directory at runtime. `-align` computes a new affine transform (channel mapper); `-transform` applies an existing one (used whenever a source-channel image needs to be warped onto the target frame).
- **Diagnostic mode:** each class has a `private final boolean diagnostic_mode = true` field controlling whether intermediate images (averages, smoothed foreground/background, per-frame background stacks) are written to disk alongside the primary outputs. There's no CLI flag for this — toggle the field directly when debugging.
- Exceptions specific to this plugin's own validation (e.g. image-size mismatches against a loaded mapping) are raised as `smFRETAnalysisException` (unchecked). Everything else is caught broadly in each `run()` and routed to `log.info(e)` + `IJ.handleException(e)` rather than propagated, since these are top-level SciJava commands.

## Build

Standard Maven build against the `pom-scijava` BOM/parent:

```
mvn clean package
```

The build resolves ImageJ/SciJava core artifacts from Maven Central plus the `scijava.public` repository declared in `pom.xml`. There is no test suite in this repository (no `src/test`), so `mvn test`/`mvn verify` runs no tests.

**Runtime dependency not in `pom.xml`:** the channel mapper requires the [TurboReg](https://imagej.net/plugins/turboreg) plugin (v2.0.1) to be present in the target Fiji/ImageJ installation's plugin path — it's invoked via reflection at runtime, not linked at compile time.

To use the plugin, install the built jar into a Fiji/ImageJ `plugins/` directory alongside TurboReg; the three commands then appear under `Plugins > smFRET`.
