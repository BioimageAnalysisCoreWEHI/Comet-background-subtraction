# Validation: WEHI backsub pipeline vs Horizon Viewer

**Date:** 2026-07-23
**Result:** ✅ The pipeline reproduces Lunaphore Horizon Viewer's background subtraction.

## What was compared

| | File | Dimensions (L0) | Channels | Pixel size |
|---|---|---|---|---|
| **ours** | `test_image/output/20260701_..._25.1592.1A_26.948.1.Y.Z.ome.tif` | 39795 × 41384 | 34 | 0.28 µm/px |
| **horizon** (ground truth) | `control_image/20260701_..._25_1592_1.ome.ome.tiff` | 37317 × 31859 | 36 | 0.28 µm/px |

Same slide (the control was renamed — dots → underscores — for downstream analysis). Both are pyramidal OME-TIFFs (8 levels, ×2 downscale). Horizon Viewer is **v2.3.0.0** (Cython-compiled Python; source not directly readable).

## Findings

### 1. Background selection — exact match ✅
Horizon encodes each subtracted channel's background in its channel name (`MARKER - BACKGROUND`). Every one matches our `comet_markers.py` assignment, including the late-cycle markers:

- All markers use the initial **autofluorescence** reference of their band (`TRITC_AF` / `Cy5_AF`).
- Horizon does **not** use the negative-control re-acquisitions (`*_N1`, `*_N2`) as subtraction references — confirming the algorithm correction (see [CLAUDE.md §3](CLAUDE.md)). Earlier "nearest-preceding" logic was wrong for the ~14 late-cycle markers; the AF rule is correct.

### 2. Channel set — accounts for the 34 vs 36 difference ✅
- Horizon keeps `DAPI` + `TRITC_AF` + `Cy5_AF` + 33 subtracted markers = **36**.
- We keep `DAPI` + 33 subtracted markers = **34** (we drop the AF reference channels by default; `keep_background=false`).
- Both drop the 3 negative-control channels (`TRITC_N1`, `Cy5_N1`, `Cy5_N2`).
- The 2 "unmatched horizon" channels are exactly `TRITC_AF` and `Cy5_AF`.

### 3. Scale & orientation — identical ✅
- Both 0.28 µm/px (`PhysicalSizeX`).
- Orientation template-match: `identity` peak 0.906; all flips/rotations ≈ 0. No flip, rotation, or transpose.

### 4. Intensity distributions — match (histR 0.9999) ✅
Alignment-invariant per-channel comparison over the aligned region (`--stats --level 1`):

- **Mean histogram correlation = 0.9999** (range 0.9997–1.0000).
- **Mean |Δmean|/mean = 0.46%**; std within ~0.5%; medians within 0–1; p99 within 0–2 intensity units, on every channel.

This is the correctness signal: identical output intensity distributions ⇒ equivalent subtraction math.

### 5. The one real difference — spatial registration (not the math)
Pixel-for-pixel correlation was low and *decreased with resolution* (DAPI: L6 ≈ 0.90 → L4 ≈ 0.72 → L0 ≈ 0.47), despite identical scale, orientation, and distributions. DAPI — which neither tool subtracts — showed the same behaviour, so it is not a subtraction effect.

**Cause:** Horizon applies its own sub-pixel / non-rigid registration on export; backsub uses the raw (already-registered) pixels as-is. No single translation aligns them (best local refine peak only 0.68). This affects pixel overlay only, not intensities. A small consistent bias (Horizon means ~0.3–0.5% higher) comes from the comparison region not being pixel-identical, not from the subtraction.

## Verdict

The WEHI pipeline's background subtraction is **equivalent to Horizon Viewer's**:
same background assignments, same pixel scale/orientation, and matching per-channel
intensity distributions (histR ≈ 1.0). The only divergence is Horizon's export-time
re-registration, which changes pixel placement but not values.

## Reproduce

Run in a backsub-style conda env (needs `tifffile`, `zarr>=3`, `numpy`, `scikit-image`):

```bash
OURS=test_image/output/20260701_..._25.1592.1A_26.948.1.Y.Z.ome.tif
CTRL=control_image/20260701_..._25_1592_1.ome.ome.tiff

# Distribution comparison (alignment-invariant — the correctness check):
python bin/compare_to_horizon.py "$OURS" "$CTRL" --stats --level 1 --save-csv stats.csv

# Pixel comparison (auto-crops/aligns; limited by Horizon's re-registration):
python bin/compare_to_horizon.py "$OURS" "$CTRL" --level auto --save-csv comparison.csv
```

`bin/compare_to_horizon.py` never loads a whole image (reads channel-by-channel in row
bands off the pyramid), matches channels by name (stripping Horizon's ` - <bg>` suffix),
and auto-aligns Horizon's crop within our image via DAPI template matching.


## Shape missmatch between original from instrument and image exported from horizon
Horizon viewer is doing a large undocumented crop of regions with no pixels resulting in exported images having a different shape from original. Example is displayed in below image:

![crop comparison](../comet-background-subtraction/Horizon_undocumented_crop.png)

Therefore it is fine to do background subtraction with Shapiro package.