#!/usr/bin/env python
'''
Module      : compare_to_horizon
Description : Pixel-wise comparison of our backsub output OME-TIFF against a
              Horizon-subtracted "ground truth" OME-TIFF, channel by channel.

              Both COMET outputs are pyramidal and huge (~100+ GB), so this NEVER
              loads a full image: it reads one channel at a time in horizontal row
              bands (via the requested pyramid level) and accumulates statistics in
              float64. Channels are matched by name (Horizon's "marker - background"
              suffix is stripped); falls back to positional matching if names don't
              line up and the channel counts are equal.

              Horizon often CROPS its export to the tissue region, so the two images
              can have different dimensions even though they are the same slide. When
              the level dimensions differ, this script auto-aligns: it template-matches
              a reference channel (default DAPI) at a coarse level to find the crop
              offset, then crops our image to Horizon's region before comparing. Pass
              --no-align to disable, or --align-channel to change the reference.

              Reported per matched channel:
                MAE       mean absolute difference (raw uint16 intensity units)
                RMSE      root-mean-square difference
                maxdiff   maximum absolute difference
                %diff     percent of pixels that differ at all
                r         Pearson correlation between the two channels
                mean_a/b  channel means (ours / horizon)

Run in the backsub conda env (needs tifffile, zarr>=3, numpy, scikit-image). The
first path is OUR PIPELINE OUTPUT (the .ome.tif in the outdir), NOT the raw input.
Set shell vars first to avoid line-continuation issues:

  OURS=test_image/output/20260701_..._Y.Z.ome.tif
  CTRL=control_image/20260701_..._25_1592_1.ome.ome.tiff
  python bin/compare_to_horizon.py "$OURS" "$CTRL" --level auto --save-csv comparison.csv

Use --level 0 for a full-resolution comparison, or --level auto (default) for a fast
downsampled check.

Copyright   : (c) WEHI SODA Hub, 2026
License     : MIT
'''

import argparse
import sys
import xml.etree.ElementTree as ET
from pathlib import Path

import numpy as np
import tifffile
import zarr

OME_NS = 'http://www.openmicroscopy.org/Schemas/OME/2016-06'


def channel_names(tif):
    """Ordered channel names from OME-XML (Channel/@Name); [] if unavailable."""
    xml = tif.ome_metadata
    if not xml:
        return []
    try:
        root = ET.fromstring(xml)
    except ET.ParseError:
        return []
    return [ch.attrib.get('Name') for ch in root.iter()
            if ch.tag.split('}', 1)[-1] == 'Channel']


def norm(name):
    """Normalise a channel name for matching: drop Horizon's ' - <bg>' suffix."""
    if name is None:
        return ''
    return name.strip().split(' - ')[0].strip().lower()


def find_channel(names, want):
    """Index of the first channel whose normalised name matches `want`, else None."""
    w = norm(want)
    for i, nm in enumerate(names):
        if norm(nm) == w:
            return i
    return None


def pick_level(series, target):
    """Pyramid-level index whose longest edge is <= target (else the smallest
    level available). Level 0 = full resolution."""
    shapes = [lvl.shape for lvl in series.levels]
    for i, shp in enumerate(shapes):
        if max(shp) <= target:
            return i
    return len(shapes) - 1


def axis_indices(axes):
    """Locate C/Y/X positions in a tifffile axes string (e.g. 'CYX')."""
    try:
        return axes.index('C'), axes.index('Y'), axes.index('X')
    except ValueError:
        raise SystemExit(f"Unsupported axis layout '{axes}' (expected C, Y and X).")


def read_band(z, ci, yi, xi, c, y0, y1, x0=0, x1=None):
    """Read channel c, rows [y0:y1), cols [x0:x1), returned as 2D (Y, X)."""
    slicer = [slice(None)] * z.ndim
    slicer[ci] = c
    slicer[yi] = slice(y0, y1)
    slicer[xi] = slice(x0, x1)
    arr = np.asarray(z[tuple(slicer)])
    return arr.T if yi > xi else arr


def open_level(series, level):
    """Return (zarr array, (ci, yi, xi)) for a pyramid level.

    For a pyramidal OME-TIFF, series.aszarr() is a multiscale Group keyed by level;
    for level 0 (the top series) an individual .aszarr() would return that Group
    rather than an array, so open the group and select the level by matching shape.
    """
    z = zarr.open(series.aszarr(), mode='r')
    if not hasattr(z, 'shape'):  # multiscale group
        target = tuple(series.levels[level].shape)
        arr = None
        for k in list(z.keys()):
            try:
                if tuple(z[k].shape) == target:
                    arr = z[k]
                    break
            except Exception:
                pass
        if arr is None:
            arr = z[str(level)]
        z = arr
    return z, axis_indices(series.levels[level].axes)


def _best_patch_center(img, box):
    """Center (y, x) of the box-sized window with the most total signal (integral
    image), so refinement lands on real structure rather than a blank region."""
    bh, bw = min(box, img.shape[0]), min(box, img.shape[1])
    ii = np.zeros((img.shape[0] + 1, img.shape[1] + 1), np.float64)
    ii[1:, 1:] = np.cumsum(np.cumsum(img.astype(np.float64), 0), 1)
    s = ii[bh:, bw:] - ii[:-bh, bw:] - ii[bh:, :-bw] + ii[:-bh, :-bw]
    y, x = np.unravel_index(int(np.argmax(s)), s.shape)
    return y + bh // 2, x + bw // 2


def align_crop(sa, sb, ref_a, ref_b, coarse_level, refine=True, patch=1024,
               search=64, min_refine=0.3):
    """Locate horizon's crop within ours using the reference channel.

    Coarse pass: template-match at `coarse_level`. Refine pass (default): match a
    full-resolution patch of horizon -- taken over its highest-signal region -- into
    a small window of ours around the coarse guess, giving an EXACT full-resolution
    integer offset (needed so nuclei line up at fine levels). If refinement is weak
    (low-texture patch), the reliable coarse offset is kept. Returns ((oy, ox)
    full-res offset, coarse_peak, refined_peak-or-None).
    """
    from skimage.feature import match_template

    za, (aci, ayi, axi) = open_level(sa, coarse_level)
    zb, (bci, byi, bxi) = open_level(sb, coarse_level)
    a2 = read_band(za, aci, ayi, axi, ref_a, 0, za.shape[ayi]).astype(np.float32)
    b2 = read_band(zb, bci, byi, bxi, ref_b, 0, zb.shape[byi]).astype(np.float32)
    if a2.shape[0] < b2.shape[0] or a2.shape[1] < b2.shape[1]:
        raise SystemExit(
            f"At align level {coarse_level} ours {a2.shape} is smaller than horizon "
            f"{b2.shape}; cannot template-match. Horizon is expected to be a crop of ours.")
    res = match_template(a2, b2)
    yx = np.unravel_index(int(np.argmax(res)), res.shape)
    scale = 2 ** coarse_level
    foy, fox = int(yx[0]) * scale, int(yx[1]) * scale
    coarse_peak = float(res.max())
    if not refine:
        return (foy, fox), coarse_peak, None

    # Refine to an exact full-res integer offset using a HIGH-SIGNAL horizon patch
    # (the center of Horizon's crop is often blank tissue -> zero-variance NCC).
    za0, (aci, ayi, axi) = open_level(sa, 0)
    zb0, (bci, byi, bxi) = open_level(sb, 0)
    Ho, Wo = za0.shape[ayi], za0.shape[axi]
    Hh, Wh = zb0.shape[byi], zb0.shape[bxi]
    ph, pw = min(patch, Hh), min(patch, Wh)
    cyc, cxc = _best_patch_center(b2, max(8, patch // scale))  # coarse coords
    cy, cx = cyc * scale, cxc * scale                          # full-res center
    py0 = min(max(0, cy - ph // 2), Hh - ph)
    px0 = min(max(0, cx - pw // 2), Wh - pw)
    tb_patch = read_band(zb0, bci, byi, bxi, ref_b, py0, py0 + ph, px0, px0 + pw).astype(np.float32)
    sy0, sx0 = max(0, foy + py0 - search), max(0, fox + px0 - search)
    sy1, sx1 = min(Ho, foy + py0 + ph + search), min(Wo, fox + px0 + pw + search)
    ta_win = read_band(za0, aci, ayi, axi, ref_a, sy0, sy1, sx0, sx1).astype(np.float32)
    if ta_win.shape[0] <= tb_patch.shape[0] or ta_win.shape[1] <= tb_patch.shape[1]:
        return (foy, fox), coarse_peak, None  # window too small to refine
    res2 = match_template(ta_win, tb_patch)
    yx2 = np.unravel_index(int(np.argmax(res2)), res2.shape)
    rpeak = float(res2.max())
    if rpeak < min_refine:
        return (foy, fox), coarse_peak, rpeak  # refine unreliable -> keep coarse
    eoy = sy0 + int(yx2[0]) - py0
    eox = sx0 + int(yx2[1]) - px0
    return (eoy, eox), coarse_peak, rpeak


def compare_channel(za, zb, ai, bi, axcols, row_chunk, oa, Hc, Wc):
    """Diff statistics for one channel pair. `oa` = (y, x) offset into ours;
    horizon is read from its origin. Comparison region is Hc x Wc."""
    (aci, ayi, axi), (bci, byi, bxi) = axcols
    oay, oax = oa

    n = 0; sum_abs = 0.0; sum_sq = 0.0; max_abs = 0; ndiff = 0
    sa = sb = saa = sbb = sab = 0.0

    for by0 in range(0, Hc, row_chunk):
        by1 = min(by0 + row_chunk, Hc)
        a = read_band(za, aci, ayi, axi, ai, oay + by0, oay + by1, oax, oax + Wc).astype(np.int64)
        b = read_band(zb, bci, byi, bxi, bi, by0, by1, 0, Wc).astype(np.int64)
        h = min(a.shape[0], b.shape[0]); w = min(a.shape[1], b.shape[1])
        a = a[:h, :w]; b = b[:h, :w]
        d = a - b; ad = np.abs(d)
        n += d.size
        sum_abs += float(ad.sum())
        sum_sq += float((d.astype(np.float64) ** 2).sum())
        m = int(ad.max()) if ad.size else 0
        if m > max_abs:
            max_abs = m
        ndiff += int(np.count_nonzero(d))
        af = a.astype(np.float64); bf = b.astype(np.float64)
        sa += af.sum(); sb += bf.sum()
        saa += (af * af).sum(); sbb += (bf * bf).sum(); sab += (af * bf).sum()

    mae = sum_abs / n if n else float('nan')
    rmse = (sum_sq / n) ** 0.5 if n else float('nan')
    pct_diff = 100.0 * ndiff / n if n else float('nan')
    mean_a = sa / n if n else float('nan')
    mean_b = sb / n if n else float('nan')
    cov = sab - sa * sb / n
    va = saa - sa * sa / n
    vb = sbb - sb * sb / n
    r = cov / (va ** 0.5 * vb ** 0.5) if va > 0 and vb > 0 else float('nan')

    return dict(pixels=n, mae=mae, rmse=rmse, maxdiff=max_abs, pct_diff=pct_diff,
                r=r, mean_a=mean_a, mean_b=mean_b)


def accumulate_hist(z, ci, yi, xi, c, oy, ox, Hc, Wc, row_chunk):
    """uint16 intensity histogram (65536 bins) over a region, read in row bands."""
    hist = np.zeros(65536, np.int64)
    for by0 in range(0, Hc, row_chunk):
        by1 = min(by0 + row_chunk, Hc)
        band = read_band(z, ci, yi, xi, c, oy + by0, oy + by1, ox, ox + Wc)
        hist += np.bincount(band.ravel().astype(np.int64), minlength=65536)[:65536]
    return hist


def hist_stats(h):
    """mean/std/median/p99/max from an intensity histogram."""
    n = int(h.sum())
    if n == 0:
        return dict(n=0, mean=float('nan'), std=float('nan'), p50=0, p99=0)
    vals = np.arange(h.size, dtype=np.float64)
    mean = float((vals * h).sum() / n)
    var = float((vals * vals * h).sum() / n - mean * mean)
    c = np.cumsum(h)
    p = lambda q: int(np.searchsorted(c, q * n))
    return dict(n=n, mean=mean, std=(var ** 0.5 if var > 0 else 0.0),
                p50=p(0.5), p99=p(0.99))


def hist_corr(ha, hb):
    """Pearson correlation between two histograms (distribution-shape similarity)."""
    a = ha.astype(np.float64) - ha.mean()
    b = hb.astype(np.float64) - hb.mean()
    d = (a * a).sum() ** 0.5 * (b * b).sum() ** 0.5
    return float((a * b).sum() / d) if d > 0 else float('nan')


def match_channels(names_a, names_b, n_a, n_b):
    """Return list of (our_idx, horizon_idx, label) and the used horizon indices."""
    pairs, used_b = [], set()
    if names_a and names_b:
        idx_b = {}
        for j, nm in enumerate(names_b):
            idx_b.setdefault(norm(nm), j)
        for i, nm in enumerate(names_a):
            j = idx_b.get(norm(nm))
            if j is not None and j not in used_b:
                used_b.add(j)
                pairs.append((i, j, nm))
    if not pairs:
        if n_a != n_b:
            raise SystemExit("No channel-name matches and channel counts differ; "
                             "cannot compare. Check the two files' channel lists.")
        print("WARNING: no channel-name matches; POSITIONAL matching.", file=sys.stderr)
        pairs = [(i, i, names_a[i] if i < len(names_a) else f'ch{i}') for i in range(n_a)]
        used_b = set(range(n_b))
    return pairs, used_b


def main():
    ap = argparse.ArgumentParser(
        description="Channel-by-channel pixel comparison of our backsub output "
                    "against a Horizon-subtracted ground-truth OME-TIFF.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    ap.add_argument('ours', type=Path, help='Our pipeline output OME-TIFF.')
    ap.add_argument('horizon', type=Path, help='Horizon-subtracted ground-truth OME-TIFF.')
    ap.add_argument('--level', default='auto',
                    help="Pyramid level to compare: integer (0 = full res) or 'auto'.")
    ap.add_argument('--target', type=int, default=4096,
                    help="For --level auto: pick the level whose longest edge is <= this.")
    ap.add_argument('--row-chunk', type=int, default=2048, help='Rows read per band.')
    ap.add_argument('--no-align', action='store_true',
                    help='Disable auto crop-alignment; require identical dimensions.')
    ap.add_argument('--align-channel', default='DAPI',
                    help='Reference channel name used to locate the crop offset.')
    ap.add_argument('--align-target', type=int, default=2048,
                    help='Coarse level (longest edge <= this) used for alignment.')
    ap.add_argument('--align-min-peak', type=float, default=0.3,
                    help='Minimum template-match correlation to accept the alignment.')
    ap.add_argument('--no-refine', action='store_true',
                    help='Skip the full-resolution offset refinement (coarse align only).')
    ap.add_argument('--channels', default=None,
                    help='Comma-separated channel names to compare (default: all matched). '
                         'Useful for a fast full-res check, e.g. --channels DAPI --level 0.')
    ap.add_argument('--stats', action='store_true',
                    help='Compare alignment-invariant intensity distributions (mean/std/'
                         'median/p99 + histogram correlation) instead of per-pixel diffs. '
                         'Use this when the two exports differ in spatial registration.')
    ap.add_argument('--save-csv', type=Path, default=None, help='Write per-channel CSV.')
    args = ap.parse_args()

    for p in (args.ours, args.horizon):
        if not p.exists():
            raise SystemExit(f"Error: file not found: {p}")

    with tifffile.TiffFile(args.ours) as ta, tifffile.TiffFile(args.horizon) as tb:
        sa, sb = ta.series[0], tb.series[0]
        names_a, names_b = channel_names(ta), channel_names(tb)

        # Choose the comparison level (same index for both; scales match at 0.28 um/px).
        if args.level == 'auto':
            lc = pick_level(sa, args.target)
        else:
            lc = int(args.level)
        lc = min(lc, len(sa.levels) - 1, len(sb.levels) - 1)

        za, axa = open_level(sa, lc)
        zb, axb = open_level(sb, lc)
        axcols = (axa, axb)
        (aci, ayi, axi), (bci, byi, bxi) = axcols
        Ho, Wo = za.shape[ayi], za.shape[axi]
        Hh, Wh = zb.shape[byi], zb.shape[bxi]
        print(f"ours   : level {lc}  (Y,X)=({Ho},{Wo})  channels={za.shape[aci]}")
        print(f"horizon: level {lc}  (Y,X)=({Hh},{Wh})  channels={zb.shape[bci]}")

        # Alignment: locate Horizon's crop within ours (unless dims already match).
        oa = (0, 0)
        if (Ho, Wo) == (Hh, Wh):
            Hc, Wc = Ho, Wo
            print("dimensions match exactly; no alignment needed.")
        elif args.no_align:
            raise SystemExit(
                f"Dimensions differ ({Ho}x{Wo} vs {Hh}x{Wh}) and --no-align is set.")
        else:
            ref_a = find_channel(names_a, args.align_channel)
            ref_b = find_channel(names_b, args.align_channel)
            if ref_a is None or ref_b is None:
                print(f"WARNING: '{args.align_channel}' not found in both "
                      f"(ours={ref_a}, horizon={ref_b}); using channel 0.", file=sys.stderr)
                ref_a = ref_a if ref_a is not None else 0
                ref_b = ref_b if ref_b is not None else 0
            la = min(pick_level(sa, args.align_target), len(sa.levels) - 1, len(sb.levels) - 1)
            (foy, fox), cpeak, rpeak = align_crop(sa, sb, ref_a, ref_b, la,
                                                  refine=not args.no_refine)
            oay = int(round(foy / 2 ** lc))
            oax = int(round(fox / 2 ** lc))
            refined_ok = rpeak is not None and rpeak >= 0.3
            rtxt = (f" refined_peak={rpeak:.4f}" if rpeak is not None else " (no refine)")
            print(f"align  : ref='{args.align_channel}'  coarse_peak={cpeak:.4f}{rtxt}  "
                  f"full-res offset (y,x)=({foy},{fox})  -> level-{lc} offset=({oay},{oax})")
            # Gate on the coarse match (reliable); refinement only adds precision.
            if cpeak < args.align_min_peak:
                raise SystemExit(
                    f"Coarse alignment correlation {cpeak:.4f} < {args.align_min_peak}; the "
                    f"crop offset is unreliable (are these really the same slide?).")
            if rpeak is not None and not refined_ok:
                print(f"WARNING: refinement weak (peak={rpeak:.4f}); using the coarse "
                      f"offset, which is only accurate to ~{2 ** la} full-res px. Expect "
                      f"residual misalignment at fine levels — compare at a coarser --level "
                      f"or pass --align-channel with a higher-signal channel.", file=sys.stderr)
            # Comparison region = Horizon's extent, clamped to ours after the offset.
            Hc = min(Hh, Ho - oay)
            Wc = min(Wh, Wo - oax)
            oa = (oay, oax)
            if Hc <= 0 or Wc <= 0:
                raise SystemExit("Computed crop falls outside ours; alignment failed.")
            print(f"       : comparing {Hc}x{Wc} region of ours (offset {oa}) vs horizon.")

        pairs, used_b = match_channels(names_a, names_b, za.shape[aci], zb.shape[bci])
        if args.channels:
            want = {norm(c) for c in args.channels.split(',')}
            pairs = [p for p in pairs if norm(p[2]) in want]
            if not pairs:
                raise SystemExit(f"None of --channels {args.channels} matched.")
        matched = {i for i, _, _ in pairs}
        unmatched_a = [n for k, n in enumerate(names_a) if k not in matched]
        unmatched_b = [n for k, n in enumerate(names_b) if k not in used_b]
        print(f"matched {len(pairs)} channels"
              + (f"; unmatched ours: {unmatched_a}" if unmatched_a else "")
              + (f"; unmatched horizon: {unmatched_b}" if unmatched_b else ""))
        print()

        if args.stats:
            sh = (f"{'channel':22s} {'mean_o':>9s} {'mean_h':>9s} {'std_o':>9s} "
                  f"{'std_h':>9s} {'p50_o':>6s} {'p50_h':>6s} {'p99_o':>6s} {'p99_h':>6s} "
                  f"{'histR':>7s}")
            print(sh)
            print('-' * len(sh))
            srows = []
            for ai, bi, label in pairs:
                ha = accumulate_hist(za, aci, ayi, axi, ai, oa[0], oa[1], Hc, Wc, args.row_chunk)
                hb = accumulate_hist(zb, bci, byi, bxi, bi, 0, 0, Hc, Wc, args.row_chunk)
                sao, sbo, hr = hist_stats(ha), hist_stats(hb), hist_corr(ha, hb)
                print(f"{str(label)[:22]:22s} {sao['mean']:9.2f} {sbo['mean']:9.2f} "
                      f"{sao['std']:9.2f} {sbo['std']:9.2f} {sao['p50']:6d} {sbo['p50']:6d} "
                      f"{sao['p99']:6d} {sbo['p99']:6d} {hr:7.4f}")
                srows.append((label, sao, sbo, hr))
            if srows:
                mr = float(np.nanmean([hr for *_, hr in srows]))
                mdev = float(np.nanmean([abs(a['mean'] - b['mean']) / max(b['mean'], 1e-9)
                                         for _, a, b, _ in srows])) * 100
                print('-' * len(sh))
                print(f"mean histogram-correlation = {mr:.4f};  mean |Δmean|/mean = {mdev:.2f}%")
                print()
                if mr >= 0.99 and mdev < 2:
                    print("VERDICT: distributions match - subtraction is equivalent to Horizon. "
                          "Remaining pixel differences are registration, not math.")
                elif mr >= 0.95:
                    print("VERDICT: distributions very close; check any channel with low histR.")
                else:
                    print("VERDICT: distributions differ - investigate channels with low histR.")
            if args.save_csv:
                import csv
                with args.save_csv.open('w', newline='') as fh:
                    w = csv.writer(fh)
                    w.writerow(['channel', 'mean_ours', 'mean_horizon', 'std_ours',
                                'std_horizon', 'p50_ours', 'p50_horizon', 'p99_ours',
                                'p99_horizon', 'hist_corr'])
                    for label, a, b, hr in srows:
                        w.writerow([label, f"{a['mean']:.4f}", f"{b['mean']:.4f}",
                                    f"{a['std']:.4f}", f"{b['std']:.4f}", a['p50'], b['p50'],
                                    a['p99'], b['p99'], f"{hr:.6f}"])
                print(f"\nWrote {args.save_csv}")
            return

        header = (f"{'channel':22s} {'MAE':>10s} {'RMSE':>10s} {'maxdiff':>8s} "
                  f"{'%diff':>7s} {'r':>8s} {'mean_ours':>10s} {'mean_horz':>10s}")
        print(header)
        print('-' * len(header))
        rows = []
        for ai, bi, label in pairs:
            s = compare_channel(za, zb, ai, bi, axcols, args.row_chunk, oa, Hc, Wc)
            print(f"{str(label)[:22]:22s} {s['mae']:10.3f} {s['rmse']:10.3f} "
                  f"{s['maxdiff']:8d} {s['pct_diff']:6.2f}% {s['r']:8.4f} "
                  f"{s['mean_a']:10.2f} {s['mean_b']:10.2f}")
            rows.append((label, s))

        if rows:
            mae_all = np.nanmean([s['mae'] for _, s in rows])
            r_all = np.nanmean([s['r'] for _, s in rows])
            max_all = max(s['maxdiff'] for _, s in rows)
            worst_r = np.nanmin([s['r'] for _, s in rows])
            print('-' * len(header))
            print(f"{'ALL (mean)':22s} {mae_all:10.3f} {'':10s} {max_all:8d} "
                  f"{'':7s} {r_all:8.4f}")
            print()
            if max_all == 0:
                print("VERDICT: identical - every matched channel is pixel-for-pixel equal.")
            elif worst_r >= 0.999 and mae_all < 1.0:
                print("VERDICT: effectively identical (r>=0.999, sub-1 MAE). Differences "
                      "are rounding/encoding-level.")
            elif worst_r >= 0.99:
                print("VERDICT: very close; inspect channels with lowest r / highest MAE. "
                      "(A small uniform MAE can also be sub-pixel crop misalignment.)")
            else:
                print("VERDICT: channels diverge (some r<0.99). Check alignment peak, "
                      "channel matching, and background assignments.")

        if args.save_csv:
            import csv
            with args.save_csv.open('w', newline='') as fh:
                w = csv.writer(fh)
                w.writerow(['channel', 'pixels', 'mae', 'rmse', 'maxdiff', 'pct_diff',
                            'pearson_r', 'mean_ours', 'mean_horizon'])
                for label, s in rows:
                    w.writerow([label, s['pixels'], f"{s['mae']:.6f}", f"{s['rmse']:.6f}",
                                s['maxdiff'], f"{s['pct_diff']:.6f}", f"{s['r']:.6f}",
                                f"{s['mean_a']:.6f}", f"{s['mean_b']:.6f}"])
            print(f"\nWrote {args.save_csv}")


if __name__ == '__main__':
    main()
