# Comet-background-subtraction — build notes

Nextflow pipeline that background-subtracts raw Lunaphore **COMET** OME-TIFFs on HPC
(parallel, one job per image) and writes pyramidal OME-TIFFs ready for QuPath — replicating
what **Horizon Viewer** does interactively in "auto" mode, but headless and batched.

Status: **scaffolded + marker detection validated on real data. Not yet run end-to-end** (Nextflow
binary not on PATH; backsub subtraction step untested here). Files written: `main.nf`,
`nextflow.config`, `nextflow_schema.json`, `bin/comet_markers.py`, `bin/environment.yml`, `README.md`,
`assets/samplesheet_example.csv`.

### Ground-truth for validation
A Horizon-subtracted version of the sample already exists in this dir (user-generated 2026-07-22):
`20260701_..._25_1592_1.ome.ome.tiff` (~117 GB). Use it to validate our backsub output pixel-wise
once the pipeline runs — do NOT delete it. (117 GB vs the 173 GB raw ≈ 39→~34 channels, consistent
with dropping the 5 AF/NC background channels; a good sanity signal.)

### Validated (stub, 2026-07-22)
`nextflow -stub-run` passes for: samplesheet input, dir-glob input (correct `sample_id` stems for
`*.ome.tif` and `*.tiff`), and both input guards (both-supplied / neither → clean error). Publishes
`<id>.ome.tif` + `markers/<id>_markers.csv` + `markers/<id>_markers_out.csv`. `comet_markers.py`
verified against the real 173 GB file (header-only, 1.4 s, exact canonical mapping).

### Remaining to run/verify
1. **First real run** — env build + backsub is untested here. `-profile conda,large` on the sample,
   then compare `<id>.ome.tif` to the Horizon ground-truth output.
2. Confirm backsub subtracts BEFORE dropping `remove=TRUE` channels (markers reference AF/NC channels
   that are themselves removed — canonical `-rr` relies on this; verify on first real run).
3. Confirm `backsub v0.5.2` pip-installs cleanly in the conda env on this cluster (first `-profile conda` run builds it).

## 9. How to run (WEHI cluster)

Nextflow is NOT on the base PATH — load it via module (needs its apptainer prereq; pins JDK 18).
`run.sh` wraps this; or manually:

```bash
module load nextflow/24.04.2                       # auto-loads apptainer, sets NXF_* + JDK18
export TMPDIR="$PWD/tmp"; mkdir -p "$TMPDIR"
nextflow run main.nf \
    --input assets/samplesheet_example.csv \
    --outdir /vast/scratch/users/$USER/comet_backsub \
    -profile conda,large
```

- For non-interactive bash, first `source /usr/local/modules/*/init/bash` before `module load`.
- Stub smoke-test (no env build, no SLURM, no pixel read): add `-stub-run` and drop the `conda` profile.
- Base process defaults are local-friendly (cpus=2); real CPU/mem come from `-profile small|medium|large`.
- `run.sh` also sets `NXF_APPTAINER_HOME_MOUNT=true` and `NXF_APPTAINER_LIBRARYDIR=/vast/projects/soda_cache/containers`.

### Locked decisions (confirmed with user 2026-07-22)
- **Input model: BOTH** — a samplesheet CSV (primary) *and* a directory-glob fallback.
- **backsub version: v0.5.2 via pip** into the conda env (packaged; identical CLI/formula to v0.4.1).
- **Drop pure background/AF channels by default** from QuPath output (via markers `remove` column);
  keep DAPI (registration/segmentation) and all subtracted markers. Flag to keep them if needed.

---

## 1. Goal (from `objectives.md`)

- User submits paths to a set of **non-subtracted** COMET images + an output dir.
- Each image is processed as an **independent SLURM job** (fan-out across nodes).
- Output = background-subtracted **pyramidal `.ome.tif`** analysable in QuPath.
- Channel → background(negative-control) mapping must be **auto-detected from the OME-XML
  metadata**, matching how Horizon Viewer's *auto* mode picks the background channel.
- Follow the **schema + profiles** style of `../export_large_annotation_regions` and the
  **Python/conda** style of `../opal_inform_stitch` (conda profile is the preferred runtime).

Sample raw image (173 GB!):
`20260701_191242_1_s9tJga_CSL_batch1_cohort_25.1592.1A_26.948.1.Y.Z.tiff`

---

## 2. The science: how background subtraction works

`backsub` (schapirolabor/background_subtraction) does **pixel-wise, exposure-scaled subtraction**:

```
Marker_corrected = Marker_raw − Background_raw × (Exposure_Marker / Exposure_Background)
```

For each marker channel it subtracts a **negative-control / autofluorescence channel** acquired
in the same spectral band, scaled by the exposure-time ratio. This is the same operation Horizon
Viewer performs; "auto mode" just means Horizon reads the channel→background mapping that the COMET
instrument wrote into the file metadata. **We read that same mapping.**

---

## 3. How channels + negative controls are detected (the OME-XML — THE crux, validated on real data)

**Validated against the real 173 GB sample** (`tiff.ome_metadata` header-only read). Real COMET
OME-XML differs from the synthetic test fixture in important ways. Real structure:

```xml
<OME xmlns="http://www.openmicroscopy.org/Schemas/OME/2016-06">
  <Image><Pixels SizeC="39" Type="uint16" DimensionOrder="XYCZT" ...>
    <Channel ID="Channel:0" Name="DAPI"     .../>   <!-- marker_name per channel -->
    <Channel ID="Channel:1" Name="TRITC_AF" .../>
    <Channel ID="Channel:3" Name="CSL3_CD68".../>
    ...
    <Plane TheC="0" ExposureTime="25"  ExposureTimeUnit="ms" .../>   <!-- exposure per channel -->
    <Plane TheC="1" ExposureTime="250" ExposureTimeUnit="ms" .../>
  </Pixels></Image>
  <StructuredAnnotations><XMLAnnotation><Value><PrivateFields ...>
    <!-- links to channel via ChannelID (NOT ID); its own ID is ChannelPriv:N -->
    <ChannelPriv ID="ChannelPriv:0" ChannelID="Channel:0" CyclePrivID="CyclePriv:0"
                 FluorescenceChannel="DAPI"  LedCurrent="20"   .../>
    <ChannelPriv ID="ChannelPriv:1" ChannelID="Channel:1" CyclePrivID="CyclePriv:0"
                 FluorescenceChannel="TRITC" LedCurrent="1700" .../>
    ...
    <!-- cycles carry the Signal/Background classification -->
    <CyclePriv ID="CyclePriv:0"  Name="Autofluorescence"          Type="AutoFluorescenceCycle" SignalType="Background"/>
    <CyclePriv ID="CyclePriv:1"  Name="Staining and Elution ..."  Type="StainingElutionCycle"  SignalType="Signal"/>
    <CyclePriv ID="CyclePriv:11" Name="Negative control"          Type="NegativeControlCycle"  SignalType="Background"/>
    ...
  </PrivateFields></Value></XMLAnnotation></StructuredAnnotations>
</OME>
```

### Key real-data facts (differ from the synthetic fixture — do not trust the fixture)
1. **`ChannelPriv` links to its channel via `@ChannelID`** (`Channel:N`); its own `@ID` is
   `ChannelPriv:N`. (Fixture wrongly used `ID="Channel:N"`.) There is **no `CycleID`** attr; it's
   `@CyclePrivID`.
2. **`FluorescenceChannel` = the spectral BAND, not a channel name** — values are `DAPI`, `TRITC`,
   `Cy5`. **No channel is literally named `TRITC` or `Cy5`.** The band's actual background/reference
   images are the autofluorescence / negative-control acquisitions: `TRITC_AF`, `Cy5_AF` (AF cycle),
   and later re-acquisitions `TRITC_N1`, `Cy5_N1`, `Cy5_N2` (NegativeControl cycles).
3. **Cycles classify Signal vs Background** via `CyclePriv/@SignalType` (`Signal` vs `Background`)
   and `@Type` (`AutoFluorescenceCycle` / `NegativeControlCycle` / `StainingElutionCycle`).

### The real "auto mode" rule (canonical — from schapirolabor's `backsub/metadata2markers.py`)
The authors of backsub (who built it for COMET) use **nearest-preceding same-band background**:
- `marker_name` ← `Channel/@Name` (deduplicate with `_1/_2` suffixes if repeated).
- `exposure`    ← `Plane/@ExposureTime` (match by `TheC`).
- `band`        ← `ChannelPriv/@FluorescenceChannel` (via `ChannelID`).
- **registration filter** (default `DAPI`): channels with `band == DAPI` are references → `background = ""`.
- A channel is a **background/reference** channel if it is in a `SignalType="Background"` cycle
  (AF or NegativeControl) — equivalently its name contains the band token (`*_AF`, `*_N1`, ...).
- For each **signal** marker (a `StainingElutionCycle` channel), `background` = **the most recent
  preceding channel of the same band that is a background/reference channel**.

Traced on the real file (verified): markers in cycles 1–10 → `TRITC_AF`/`Cy5_AF`; cycles 12–19 →
`TRITC_N1`/`Cy5_N1`; cycles 21–22 → `Cy5_N2` (TRITC in those late cycles → nearest prior `TRITC_N1`).
This is what Horizon "auto" effectively does: pick the same-band blank acquired closest before the marker.

> ⚠️ **`../sp_segment/bin/extract_markers.py` is INSUFFICIENT for this real data.** It assumes
> `FluorescenceChannel` equals a real channel name and would emit `background=TRITC`/`Cy5` (nonexistent
> channels) → backsub fails. It also matched `ChannelPriv/@ID` (never `ChannelID`), so on real files it
> maps nothing. **We implement the canonical nearest-preceding algorithm instead**, in our own tested
> `bin/comet_markers.py`. (`metadata2markers.py`'s own `main()` is also out of sync with its argparse /
> function signatures and isn't a console entry point — port the `assign_background` logic, don't call it.)

Header-only read only — the file is 173 GB; never `tiff.asarray()`, only `tiff.ome_metadata`.

### backsub markers CSV format (what `bin/comet_markers.py` must emit)
Columns: `marker_name, background, exposure` (+ optional `remove` = `TRUE` to drop a channel from
output). `background` must equal another row's `marker_name` exactly, or be empty. Real-file excerpt
(background channels dropped by default per locked decision):

```csv
marker_name,background,exposure,remove
DAPI,,25.0,
TRITC_AF,,250.0,TRUE
Cy5_AF,,400.0,TRUE
CSL3_CD68,TRITC_AF,250.0,
CSL_CD204,Cy5_AF,400.0,
...
CSL_FN1,TRITC_N1,250.0,       # cycle 12 → uses the re-acquired TRITC_N1 blank
CSL_CD90,Cy5_N1,400.0,
...
CSL3_GPA,Cy5_N2,400.0,        # cycle 21 → uses Cy5_N2
```

---

## 4. The `backsub` tool (schapirolabor/background_subtraction)

- Repo: https://github.com/schapirolabor/background_subtraction  — container
  `ghcr.io/schapirolabor/background_subtraction:v0.4.1` (used by sp_segment).
- CLI (stable across v0.4.x–v0.5.x):
  ```
  backsub  -r/--input <in.tiff>  -m/--markers <markers.csv> \
           -o/--output <out.ome.tif>  -mo/--marker-output <out_markers.csv> \
           [-mpp pixel_size] [-ts tile_size(mult 16, def 256)] [-dsf downscale(def 2)] \
           [-sr/--save-ram] [-comp lzw|zlib|deflate|none (def zlib)]
  ```
  (v0.4.1 container invokes `python3 /background_subtraction/background_sub.py ...`; the packaged
  v0.5.x exposes the `backsub` console command — same args.)
- Writes a **pyramidal OME-TIFF**; channels marked `remove` are dropped from output.
- Deps (from pyproject): dask, dask-image, imagecodecs, numpy, ome_types, pandas, psutil,
  scikit-image, tifffile, zarr>=3, (loguru). Python ≥3.11.

### ⚠️ Packaging / version note (drives the conda decision)
- **v0.4.1 has NO `pyproject.toml`** → not cleanly pip-installable (that's why sp_segment's
  backsub module *hard-errors under `-profile conda/mamba`* and forces a container).
- Packaging + `[project.scripts] backsub = backsub.background_sub:main` **starts at v0.5.x**
  (latest **v0.5.2**). So for a **conda profile we pip-install v0.5.2 from git**, giving the
  `backsub` command inside our own env. CLI + formula are unchanged from v0.4.1.
- Newest `backsub` also bundles `backsub/metadata2markers.py` (its own markers generator using
  `ome_types` + a `registration_filter` + substring matching). We prefer WEHI's `extract_markers.py`
  because it decodes `FluorescenceChannel` directly (closer to Horizon auto mode, simpler deps).

---

## 5. Planned architecture

Two processes per image (nf-core idiom; keeps the markers CSV inspectable/publishable for QC):

```
 input images ──► EXTRACT_MARKERS ──► markers.csv ──► BACKSUB ──► <name>.ome.tif (+ out_markers.csv)
   (parallel, one (SLURM) job per image)
```

- **EXTRACT_MARKERS** (`label process_low`): `bin/comet_markers.py <tiff> [--remove-marker X]
  [--keep-background] [--registration-filter DAPI]` → `<id>_markers.csv`. Deps: tifffile, pandas
  (stdlib ElementTree for XML). Header-only OME-XML read. Implements the §3 canonical algorithm.
- **BACKSUB** (`label process_backsub`, big RAM/time): `backsub -r <tiff> -m <markers.csv>
  -o <id>.ome.tif -mo <id>_markers_out.csv [-sr] [-comp lzw] ...`.

Input model (LOCKED: both):
- **(B) samplesheet CSV** (primary): `sample_id,image_path[,remove]` for provenance + per-image
  remove lists. `--input` / `--samplesheet`.
- **(A) directory glob** (fallback convenience): `--input_dir` → glob `*.tif/*.tiff/*.ome.tif`, fan
  out one job per file (sample_id = filename stem), like opal's `base_dir`.
- Exactly one of the two must be supplied.

Publish: `${outdir}/<id>.ome.tif`, `${outdir}/markers/<id>_markers.csv`, plus
`pipeline_info/` timeline/report/trace/dag (copy opal's config block).

### Memory reality (173 GB inputs)
backsub loads/pyramids large arrays — profiles must be generous and use `-sr/--save-ram`.
Size `process_backsub` like opal's `large` (≈960 GB) as the default target; expose `--save_ram`,
`--compression`, `--tile_size` as params. Use `task.attempt` memory scaling + retry.

---

## 6. Repo layout to create (mirror opal_inform_stitch)

```
Comet-background-subtraction/
  main.nf                      # workflow: discover images → EXTRACT_MARKERS → BACKSUB
  nextflow.config              # params, profiles, env, shell, timeline/report/trace/dag, manifest
  nextflow_schema.json         # nf-core-style schema (input, outdir, remove markers, backsub opts)
  bin/
    comet_markers.py           # canonical §3 algorithm (nearest-preceding same-band background)
    environment.yml            # conda env (see below)
  README.md
  CLAUDE.md                    # (this file)
```

### `bin/environment.yml` (conda profile)
```yaml
name: comet_backsub
channels: [conda-forge, bioconda]
dependencies:
  - python=3.12
  - tifffile
  - pandas
  - typer
  - numpy
  - scikit-image
  - dask
  - dask-image
  - ome-types
  - zarr>=3
  - imagecodecs
  - psutil
  - loguru
  - pip
  - pip:
    - "backsub @ git+https://github.com/schapirolabor/background_subtraction.git@v0.5.2"
```

### Profiles (copy opal's set)
`conda` (primary), `mamba`, `singularity`/`apptainer` (could instead use the ghcr v0.4.1 container
for BACKSUB — but then invocation path differs; keep conda as the single-env happy path), plus
resource profiles `small` / `medium` / `large` (SLURM, `queue='regular'`). Reuse opal's `env{}` block
(`PYTHONNOUSERSITE=1`, `TMPDIR`) and `process.shell` bash-strict settings.

---

## 7. Open questions (remaining) + resolved

Resolved (see Locked decisions at top): input model = both; backsub = v0.5.2 pip; drop background
channels by default. Real metadata validated (§3).

Remaining to confirm:
1. **compression**: QuPath-friendly `lzw` vs backsub default `zlib`? (export_large uses LZW.)
   Default proposal: `lzw`, param-exposed.
2. **Multiple DAPI**: this sample has one DAPI; if other slides re-acquire DAPI per cycle, decide
   whether to drop extra DAPI (canonical `--remove_dapi` keeps only the first). Expose a flag.
3. **Which same-band background is "correct"** for late TRITC markers when only `TRITC_N1` exists
   (no `TRITC_N2`): canonical picks nearest *preceding* → `TRITC_N1`. Confirm that matches Horizon.
4. **Run environment**: Nextflow binary not yet on PATH — install/locate before first real run.

---

## 8. Key reference paths (all local)

- `../export_large_annotation_regions/`  — schema.json + profiles style to follow.
- `../opal_inform_stitch/`                — Python+conda pipeline template (main.nf, config, env.yml).
- `../sp_segment/bin/extract_markers.py`  — **the auto-detection logic to reuse.**
- `../sp_segment/bin/combine_channels.py` — robust multi-format OME/ImageJ/MIBI metadata parsing patterns.
- `../sp_segment/modules/nf-core/backsub/main.nf`        — how backsub is invoked (container path).
- `../sp_segment/modules/local/extractmarkers/main.nf`   — conda+container module wiring for markers.
- `../sp_segment/tests/data/comet/make_comet_test_data.py` — synthetic COMET OME-XML generator (test fixture).

Tools present: Java 25 at `/vast/scratch/users/mckay.m/java/jdk-25.0.2+10`. **Nextflow not on PATH**
— locate/install before running. `tifffile` available in env
`/vast/scratch/users/mckay.m/envs/cellmeasurement-test` (handy for header inspection).
</content>
