# Comet-background-subtraction

Nextflow pipeline that background-subtracts raw **Lunaphore COMET** OME-TIFFs headlessly on HPC —
reproducing what **Horizon Viewer** does interactively in *auto* mode — and writes pyramidal
OME-TIFFs ready for **QuPath**. Each image is processed as an independent SLURM job, so a whole
batch is subtracted in parallel.

## What it does

For each input image:

1. **`EXTRACT_MARKERS`** reads the COMET **OME-XML metadata** (header only — the pixels are never
   loaded) and auto-derives a [backsub](https://github.com/schapirolabor/background_subtraction)
   markers CSV. It reproduces Horizon's *auto* background detection: each signal marker is paired
   with the **most recent preceding autofluorescence / negative-control channel acquired in the same
   spectral band** (`FluorescenceChannel`). Registration channels (default `DAPI`) are never
   subtracted, and pure background channels are dropped from the output by default.
2. **`BACKSUB`** runs background subtraction
   (`Marker_corrected = Marker_raw − Background × (Exposure_Marker / Exposure_Background)`) and
   writes a pyramidal OME-TIFF plus the resolved markers CSV.

### How channels + negative controls are detected

COMET writes, per channel, into the OME-XML:

- `Channel/@Name` → marker name
- `Plane/@ExposureTime` → exposure (matched by `TheC`)
- private `ChannelPriv/@FluorescenceChannel` (linked via `@ChannelID`) → the spectral **band**
  (`DAPI`/`TRITC`/`Cy5`)
- `CyclePriv/@SignalType` (`Signal` vs `Background`) → which cycles are autofluorescence /
  negative-control acquisitions

The pipeline uses this to pair each marker with the correct band-matched blank. See
[CLAUDE.md](CLAUDE.md) for the full metadata analysis and worked example on a real 39-channel slide.

## Requirements

- [Nextflow](https://www.nextflow.io/) ≥ 24.04.2
- conda or mamba (the `conda` profile builds an environment that includes `backsub` v0.5.2)

## Inputs

Provide **exactly one** of:

**A samplesheet CSV** (`--input`) — recommended:

```csv
sample_id,image_path,remove
slideA,/vast/.../slideA_raw.ome.tiff,
slideB,/vast/.../slideB_raw.ome.tiff,CSL_panCK;CSL_FAP
```

- `sample_id` — output basename (`<sample_id>.ome.tif`).
- `image_path` — absolute path to the raw (non-subtracted) COMET image.
- `remove` *(optional)* — `;` or `|` separated marker names to drop from that image.

**Or a directory** (`--input_dir`) of `*.tif` / `*.tiff` / `*.ome.tif`; one job per file, with
`sample_id` taken from the filename stem.

## Usage

```bash
nextflow run main.nf \
    --input samplesheet.csv \
    --outdir /vast/scratch/users/$USER/comet_backsub \
    -profile conda,large
```

or, over a directory:

```bash
nextflow run main.nf \
    --input_dir /path/to/raw_comet_images \
    --outdir /path/to/output \
    -profile conda,medium
```

## Parameters

| Parameter               | Default  | Description                                                                 |
|-------------------------|----------|-----------------------------------------------------------------------------|
| `--input`               | —        | Samplesheet CSV (`sample_id,image_path[,remove]`).                          |
| `--input_dir`           | —        | Directory of raw COMET images (alternative to `--input`).                   |
| `--outdir`              | —        | **Required.** Output directory.                                            |
| `--registration_filter` | `DAPI`   | Band used for registration; never subtracted.                              |
| `--keep_background`     | `false`  | Keep AF / negative-control channels in the output.                         |
| `--remove_extra_dapi`   | `false`  | Drop every registration (DAPI) channel except the first.                   |
| `--save_ram`            | `true`   | Pass `-sr` to backsub (~50% less RAM; recommended for large slides).       |
| `--compression`         | `lzw`    | Output compression: `lzw` / `zlib` / `deflate` / `none`.                   |
| `--tile_size`           | `256`    | Pyramid tile size (multiple of 16).                                        |
| `--downscale_factor`    | `2`      | Pyramid downscale factor (non-pyramidal inputs only).                      |
| `--pixel_size`          | —        | Microns/pixel; taken from metadata if unset.                               |
| `--publish_dir_mode`    | `copy`   | Nextflow `publishDir` mode.                                                |

## Resource profiles

Combine an environment profile (`conda` / `mamba` / `singularity` / `apptainer`) with a resource
profile:

| Profile  | `process_backsub` memory | Target nodes    | Use case                       |
|----------|--------------------------|-----------------|--------------------------------|
| `small`  | 128 GB                   | sml / med       | Small slides / test runs       |
| `medium` | 480 GB                   | med / il        | Medium slides                  |
| `large`  | 960 GB                   | lrg             | Full-size COMET exports (100s of GB) |

```bash
-profile conda,large
```

Memory scales with `task.attempt`, so a job that hits an OOM retries once at double the memory.

## Output

```
outdir/
    <sample_id>.ome.tif                 # background-subtracted pyramidal OME-TIFF (open in QuPath)
    markers/
        <sample_id>_markers.csv         # auto-detected markers (input to backsub)
        <sample_id>_markers_out.csv     # markers matching the output channels
    pipeline_info/
        execution_timeline_*.html
        execution_report_*.html
        execution_trace_*.txt
        pipeline_dag_*.html
```

## Notes

- `bin/comet_markers.py` can be run standalone to preview the detected markers before a full run:
  ```bash
  python bin/comet_markers.py /path/to/raw_comet.ome.tiff        # prints CSV to stdout
  python bin/comet_markers.py /path/to/raw_comet.ome.tiff --keep-background
  ```
- `backsub` has no conda/bioconda package; the `conda` profile pip-installs the packaged **v0.5.2**
  release (identical CLI and subtraction formula to the `ghcr.io/schapirolabor/background_subtraction:v0.4.1`
  container used elsewhere).
