#!/usr/bin/env python
'''
Module      : comet_markers
Description : Generate a background_subtraction (backsub) compatible markers CSV
              from Lunaphore COMET OME-TIFF metadata.

              Replicates Horizon Viewer "auto" background detection: every signal
              marker is paired with the most recent PRECEDING autofluorescence /
              negative-control channel acquired in the SAME spectral band
              (OME-XML ChannelPriv/@FluorescenceChannel). Reference / registration
              channels (default DAPI) are never subtracted.

              The pairing algorithm is a faithful port of the canonical COMET
              logic in schapirolabor/background_subtraction
              (backsub/metadata2markers.py -> assign_background), adapted to read
              the real COMET private-field schema:
                * ChannelPriv is linked to its Channel via @ChannelID
                  (its own @ID is "ChannelPriv:N", NOT "Channel:N");
                * FluorescenceChannel holds the spectral BAND (DAPI/TRITC/Cy5),
                  not a channel name;
                * CyclePriv/@SignalType marks whole cycles as "Signal" or
                  "Background" (AutoFluorescenceCycle / NegativeControlCycle).

Copyright   : (c) WEHI SODA Hub, 2026
License     : MIT
Portability : POSIX

Output columns (backsub-compatible): marker_name,background,exposure,remove
Writes to stdout by default.
'''

import argparse
import sys
import xml.etree.ElementTree as ET
from pathlib import Path

import pandas as pd
from tifffile import TiffFile

OME_NS = 'http://www.openmicroscopy.org/Schemas/OME/2016-06'
NS = {'ome': OME_NS}


def _localname(tag):
    """Strip the XML namespace from a tag, e.g. '{...}Channel' -> 'Channel'."""
    return tag.split('}', 1)[-1]


def _make_unique(names):
    """
    Ensure marker names are unique; duplicates get an index suffix (_1, _2, ...).
    Mirrors backsub.metadata2markers.make_marker_names_unique_list so that a
    background reference always resolves to exactly one row.
    """
    seen = {}
    out = []
    for name in names:
        if name in seen:
            seen[name] += 1
            out.append(f"{name}_{seen[name]}")
        else:
            seen[name] = 0
            out.append(name)
    return out


def parse_comet_metadata(ome_xml):
    """
    Parse COMET OME-XML into an ordered list of per-channel dicts with keys:
      index, marker_name, band, exposure, cycle_id, signal_type, is_background

    band          : FluorescenceChannel (spectral band, e.g. DAPI/TRITC/Cy5) or None
    signal_type   : owning cycle's SignalType ("Signal"/"Background") or None
    is_background : True if the channel is an AF / negative-control acquisition
    """
    root = ET.fromstring(ome_xml)

    # --- Channels (document order defines channel index) ---
    channels = root.findall('.//ome:Channel', NS)
    raw_names = [ch.attrib.get('Name') or ch.attrib.get('ID', f'Channel_{i}')
                 for i, ch in enumerate(channels)]
    names = _make_unique([str(n) for n in raw_names])
    channel_ids = [ch.attrib.get('ID') for ch in channels]

    # --- Exposure times, matched by Plane/@TheC (never assume plane order) ---
    exposure_by_c = {}
    for pl in root.findall('.//ome:Plane', NS):
        the_c = pl.attrib.get('TheC')
        exp = pl.attrib.get('ExposureTime')
        if the_c is not None and exp is not None:
            try:
                exposure_by_c[int(the_c)] = float(exp)
            except (TypeError, ValueError):
                pass

    # --- Private fields: iterate namespace-agnostically (they inherit the OME ns) ---
    # ChannelPriv links to its channel via @ChannelID; carries band + owning cycle.
    channelpriv_by_chid = {}
    cyclepriv = {}
    for el in root.iter():
        name = _localname(el.tag)
        if name == 'ChannelPriv':
            chid = el.attrib.get('ChannelID') or el.attrib.get('ID')
            channelpriv_by_chid[chid] = el.attrib
        elif name == 'CyclePriv':
            cyclepriv[el.attrib.get('ID')] = el.attrib

    records = []
    for i, (cid, marker) in enumerate(zip(channel_ids, names)):
        priv = channelpriv_by_chid.get(cid, {})
        band = priv.get('FluorescenceChannel')
        cycle_id = priv.get('CyclePrivID')
        cyc = cyclepriv.get(cycle_id, {})
        signal_type = cyc.get('SignalType')
        records.append({
            'index': i,
            'marker_name': marker,
            'band': band,
            'exposure': exposure_by_c.get(i),
            'cycle_id': cycle_id,
            'signal_type': signal_type,
            'is_background': False,  # filled in by classify_channels
        })
    return records


def classify_channels(records, registration_filter):
    """
    Mark which channels are background / reference (autofluorescence or
    negative-control) acquisitions and therefore not subtracted themselves.

    Primary signal: the owning cycle's SignalType == "Background".
    Fallback (when cycle metadata is absent): the channel name contains its own
    band token as a substring (e.g. TRITC_AF, Cy5_N1) -- the same heuristic the
    canonical backsub tool uses.
    """
    for r in records:
        band = r['band']
        if band is None:
            r['is_background'] = False
            continue
        if r['signal_type'] is not None:
            r['is_background'] = (r['signal_type'].lower() == 'background'
                                  and band != registration_filter)
        else:
            # Fallback: name contains the band token but is not a registration channel
            r['is_background'] = (band != registration_filter
                                  and band in r['marker_name'])
    return records


def assign_backgrounds(records, registration_filter):
    """
    For each signal marker, assign background = the most recent PRECEDING
    background/reference channel of the same band. Reference (registration) and
    background channels themselves get an empty background.

    Returns the number of signal markers left without a background (warned).
    """
    unmatched = 0
    for r in records:
        band = r['band']
        # Reference/registration channels and the background channels themselves
        # are never subtracted.
        if band is None or band == registration_filter or r['is_background']:
            r['background'] = ''
            continue
        # Signal marker: search preceding channels (nearest first) for a same-band
        # background/reference acquisition.
        bg = ''
        for prev in reversed(records[:r['index']]):
            if prev['band'] == band and prev['is_background']:
                bg = prev['marker_name']
                break
        if not bg:
            unmatched += 1
            sys.stderr.write(
                f"Warning: no preceding '{band}' background channel found for "
                f"marker '{r['marker_name']}'; leaving background empty.\n"
            )
        r['background'] = bg
    return unmatched


def build_dataframe(records, registration_filter, drop_background,
                    remove_markers, remove_extra_dapi):
    """Assemble the backsub markers DataFrame with the `remove` column."""
    remove_set = set(remove_markers or [])
    seen_registration = False
    rows = []
    for r in records:
        remove = ''
        band = r['band']
        is_registration = (band is not None and band == registration_filter)

        if r['marker_name'] in remove_set:
            remove = 'TRUE'
        elif drop_background and r['is_background']:
            # Drop pure autofluorescence / negative-control channels from output.
            remove = 'TRUE'
        elif is_registration:
            if seen_registration and remove_extra_dapi:
                remove = 'TRUE'  # keep only the first registration (DAPI) channel
            seen_registration = True

        exposure = r['exposure'] if r['exposure'] is not None else 1.0
        rows.append({
            'marker_name': r['marker_name'],
            'background': r.get('background', ''),
            'exposure': float(exposure),
            'remove': remove,
        })
    return pd.DataFrame(rows, columns=['marker_name', 'background', 'exposure', 'remove'])


def imagej_fallback(tiff, registration_filter):
    """
    When there is no OME-XML, fall back to ImageJ channel labels: no background
    detection is possible, so emit blank backgrounds and default exposure.
    """
    imagej = getattr(tiff, 'imagej_metadata', None)
    if not isinstance(imagej, dict):
        return None
    labels = imagej.get('Labels')
    if not labels:
        return None
    sys.stderr.write(
        "Warning: no OME-XML found; using ImageJ channel labels with blank "
        "backgrounds and exposure=1.0 (no background subtraction will occur).\n"
    )
    names = _make_unique([str(x) for x in labels])
    return pd.DataFrame(
        [{'marker_name': n, 'background': '', 'exposure': 1.0, 'remove': ''}
         for n in names],
        columns=['marker_name', 'background', 'exposure', 'remove'],
    )


def main():
    ap = argparse.ArgumentParser(
        description=("Generate a backsub-compatible markers CSV from a Lunaphore "
                     "COMET OME-TIFF (Horizon auto-mode background detection)."),
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    ap.add_argument('tiff', type=Path, help='Path to the COMET OME-TIFF (header only is read).')
    ap.add_argument('-o', '--output', type=Path, default=None,
                    help='Output CSV path (default: stdout).')
    ap.add_argument('-rf', '--registration-filter', default='DAPI',
                    help='FluorescenceChannel/band used for registration; these '
                         'channels are references and never subtracted.')
    ap.add_argument('-r', '--remove-marker', action='append', default=[],
                    dest='remove_marker', metavar='NAME',
                    help='Marker name to flag remove=TRUE (repeatable).')
    ap.add_argument('--keep-background', action='store_true',
                    help='Keep autofluorescence / negative-control channels in the '
                         'output (default: drop them via remove=TRUE).')
    ap.add_argument('--remove-extra-dapi', action='store_true',
                    help='Flag every registration (DAPI) channel except the first '
                         'with remove=TRUE.')
    args = ap.parse_args()

    with TiffFile(args.tiff) as tiff:
        ome_xml = tiff.ome_metadata
        if not ome_xml:
            desc = tiff.pages[0].description
            if isinstance(desc, str) and desc.lstrip().startswith('<'):
                ome_xml = desc

        df = None
        if ome_xml:
            try:
                records = parse_comet_metadata(ome_xml)
                records = classify_channels(records, args.registration_filter)
                assign_backgrounds(records, args.registration_filter)
                df = build_dataframe(
                    records,
                    registration_filter=args.registration_filter,
                    drop_background=not args.keep_background,
                    remove_markers=args.remove_marker,
                    remove_extra_dapi=args.remove_extra_dapi,
                )
            except ET.ParseError as exc:
                sys.stderr.write(f"Warning: could not parse OME-XML ({exc}).\n")
                df = None

        if df is None:
            df = imagej_fallback(tiff, args.registration_filter)

    if df is None or df.empty:
        raise SystemExit(
            "Error: could not extract channel metadata from OME-XML or ImageJ labels."
        )

    # Validate requested removals exist
    known = set(df['marker_name'])
    for name in args.remove_marker:
        if name not in known:
            raise SystemExit(f"Error: --remove-marker '{name}' not found in channels.")

    # Sanity check: every non-empty background must resolve to a real marker name.
    for bg in df.loc[df['background'] != '', 'background']:
        if bg not in known:
            raise SystemExit(
                f"Error: background '{bg}' does not match any marker_name; "
                "metadata may be non-standard."
            )

    out = args.output.open('w') if args.output else sys.stdout
    try:
        df.to_csv(out, index=False)
    finally:
        if args.output:
            out.close()


if __name__ == '__main__':
    main()
