nextflow.enable.dsl = 2

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    BioimageAnalysisCoreWEHI/Comet-background-subtraction
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Background-subtract raw Lunaphore COMET OME-TIFFs (Horizon "auto" mode, headless)
    and write pyramidal OME-TIFFs ready for QuPath. One SLURM job per image.

    Per image:
      EXTRACT_MARKERS  reads the COMET OME-XML and auto-derives the backsub markers CSV
                       (each signal marker -> nearest preceding same-band AF/negative
                       -control channel), then
      BACKSUB          runs schapirolabor/background_subtraction to produce the
                       corrected pyramidal OME-TIFF.
----------------------------------------------------------------------------------------
*/

/*
 * Derive a backsub-compatible markers CSV from the COMET OME-XML metadata.
 * Header-only read of the (potentially very large) input image.
 */
process EXTRACT_MARKERS {
    tag "${sample_id}"
    label 'process_low'

    conda "${projectDir}/bin/environment.yml"

    publishDir "${params.outdir}/markers", mode: params.publish_dir_mode, pattern: "*_markers.csv"

    input:
    tuple val(sample_id), path(image), val(remove_markers)

    output:
    tuple val(sample_id), path(image), path("${sample_id}_markers.csv"), path("${sample_id}_mpp.txt"), emit: markers

    script:
    def reg     = params.registration_filter ? "--registration-filter '${params.registration_filter}'" : ""
    def keepbg  = params.keep_background      ? "--keep-background"    : ""
    def rmdapi  = params.remove_extra_dapi    ? "--remove-extra-dapi"  : ""
    def nearest = params.use_nearest_background ? "--use-nearest-background" : ""
    def rm_args = (remove_markers ?: []).collect { "--remove-marker '${it}'" }.join(' ')
    """
    python ${projectDir}/bin/comet_markers.py \\
        "${image}" \\
        ${reg} ${keepbg} ${rmdapi} ${nearest} ${rm_args} \\
        --pixel-size-out "${sample_id}_mpp.txt" \\
        -o "${sample_id}_markers.csv"
    """

    stub:
    """
    printf 'marker_name,background,exposure,remove\\nDAPI,,25.0,\\n' > ${sample_id}_markers.csv
    printf '0.28' > ${sample_id}_mpp.txt
    """
}

/*
 * Run background subtraction, producing a pyramidal OME-TIFF for QuPath.
 */
process BACKSUB {
    tag "${sample_id}"
    label 'process_backsub'

    conda "${projectDir}/bin/environment.yml"

    publishDir "${params.outdir}",         mode: params.publish_dir_mode, pattern: "*.ome.tif"
    publishDir "${params.outdir}/markers", mode: params.publish_dir_mode, pattern: "*_markers_out.csv"

    input:
    tuple val(sample_id), path(image), path(markers), path(mpp_file)

    output:
    tuple val(sample_id), path("${sample_id}.ome.tif"), emit: image
    tuple val(sample_id), path("${sample_id}_markers_out.csv"), emit: markers
    path "versions.yml", emit: versions

    script:
    def save_ram  = params.save_ram        ? "-sr"                          : ""
    def comp      = params.compression     ? "-comp ${params.compression}"  : ""
    def tile      = params.tile_size       ? "-ts ${params.tile_size}"      : ""
    def dsf       = params.downscale_factor ? "-dsf ${params.downscale_factor}" : ""
    // Explicit --pixel_size param overrides the value detected from metadata.
    def mppParam  = params.pixel_size ? "${params.pixel_size}" : ""
    if ("${image}" == "${sample_id}.ome.tif")
        error "Input and output names collide for '${sample_id}'; rename the input or set a different sample_id."
    """
    # Pixel size: explicit --pixel_size wins; otherwise use the value comet_markers.py
    # detected from the OME metadata (written to ${mpp_file}). Passed to backsub as -mpp
    # so the micron scale is explicit and logged rather than silently re-read.
    MPP="${mppParam}"
    if [ -z "\$MPP" ] && [ -s "${mpp_file}" ]; then MPP=\$(tr -d '[:space:]' < "${mpp_file}"); fi
    if [ -n "\$MPP" ]; then
        MPP_ARG="-mpp \$MPP"
        echo "BACKSUB using pixel size: -mpp \$MPP (micrometres)"
    else
        MPP_ARG=""
        echo "BACKSUB: no pixel size available; backsub will read metadata / fall back to 1 pixel/unit"
    fi

    backsub \\
        -r "${image}" \\
        -m "${markers}" \\
        -o "${sample_id}.ome.tif" \\
        -mo "${sample_id}_markers_out.csv" \\
        ${save_ram} ${comp} ${tile} ${dsf} \$MPP_ARG

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        backsub: \$(backsub --version 2>&1 | tr -dc '0-9.')
    END_VERSIONS
    """

    stub:
    """
    touch ${sample_id}.ome.tif
    cp ${markers} ${sample_id}_markers_out.csv
    echo '"${task.process}": {backsub: stub}' > versions.yml
    """
}

workflow {
    if (!params.outdir) {
        error "Missing required parameter: --outdir"
    }

    def hasSheet = params.input      as boolean
    def hasDir   = params.input_dir  as boolean
    if (hasSheet == hasDir) {
        error "Provide exactly one of --input (samplesheet CSV) or --input_dir (directory of images)."
    }

    if (hasSheet) {
        // Samplesheet: sample_id,image_path[,remove]   (remove = ';' or '|' separated marker names)
        images = Channel.fromPath(params.input, checkIfExists: true)
            .splitCsv(header: true)
            .map { row ->
                def sid = row.sample_id?.trim()
                def ipath = row.image_path?.trim()
                if (!sid)   error "Samplesheet row missing 'sample_id': ${row}"
                if (!ipath) error "Samplesheet row missing 'image_path': ${row}"
                def img = file(ipath, checkIfExists: true)
                def rm = row.remove?.trim()
                    ? row.remove.trim().split(/[;|]/).collect { it.trim() }.findAll { it }
                    : []
                tuple(sid, img, rm)
            }
    } else {
        // Directory glob convenience: one job per *.tif / *.tiff / *.ome.tif
        images = Channel.fromPath("${params.input_dir}/*.{tif,tiff}", checkIfExists: true)
            .map { img ->
                def id = img.name.replaceAll(/\.ome\.tif(f)?$/, '').replaceAll(/\.tif(f)?$/, '')
                tuple(id, img, [])
            }
    }

    EXTRACT_MARKERS(images)
    BACKSUB(EXTRACT_MARKERS.out.markers)
}
