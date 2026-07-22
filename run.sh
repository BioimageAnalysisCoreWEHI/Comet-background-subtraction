#!/usr/bin/env bash
# Launcher for the Comet-background-subtraction pipeline on the WEHI cluster.
# Usage:
#   ./run.sh --input assets/samplesheet_example.csv --outdir /path/to/out -profile conda,large
#   ./run.sh --input_dir /path/to/images          --outdir /path/to/out -profile conda,medium
# Any extra args are passed straight through to `nextflow run`.
set -euo pipefail

module load nextflow/24.04.2

# Apptainer/Singularity settings (harmless under the conda profile)
export NXF_APPTAINER_HOME_MOUNT=true
export NXF_APPTAINER_LIBRARYDIR=/vast/projects/soda_cache/containers

# Keep temp local to the launch dir
export TMPDIR="$PWD/tmp"
mkdir -p "$TMPDIR"
echo "[INFO] Set TMPDIR to: $TMPDIR"

nextflow run "$(dirname "$0")/main.nf" "$@"
