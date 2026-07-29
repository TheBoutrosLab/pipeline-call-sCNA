#!/usr/bin/env nextflow

log.info """\
------------------------------------
        B P G - B A F / Log R
------------------------------------
Docker Images:
- docker_image_bpg: ${params.docker_image_bpg}
"""

include { generate_standard_filename; sanitize_string } from '../external/pipeline-Nextflow-module/modules/common/generate_standardized_filename/main.nf'

/**
* Plot Battenberg BAF and LogR values with BoutrosLabPlottingGeneral.
*/
process plot_BAFLogR_BPG {
    container params.docker_image_bpg

    tag "${sample_id}-${per_chrom_mode}"

    publishDir "${META.workflow_output_dir}/QC",
        pattern: "*.png",
        mode: "copy",
        saveAs: {
            "${output_filename}_${sanitize_string(file(it).getName().replace("${sample_id}_", ""))}"
        }

    ext log_dir_suffix: { "-${sample_id}-${per_chrom_mode}" }

    input:
        val(META)
        val(position_scale)
        val(tumor_id)
        tuple val(sample_id), path(baf_file), path(logr_file)
        each per_chrom_mode
        path(script_dir)
        val(plot_resolution)
        path(reference_dict)

    output:
        path "${sample_id}_*.png", emit: baf_logr_plots

    script:
    output_filename = generate_standard_filename(
        "BPG-${params.bpg_version}",
        params.dataset_id,
        sample_id,
        [additional_tools:["Battenberg-${params.battenberg_version}"]]
    )

    per_chrom_mode_arg = "--per-chrom-mode"
    reference_dict_arg = "-d ${reference_dict}"
    position_scale_arg = ""

    if ( per_chrom_mode == "disabled" ) {
        per_chrom_mode_arg = ""
        position_scale_arg = "-p ${position_scale}"
    }

    """
    set -euo pipefail

    Rscript ${script_dir}/BPG_BAF_LogR.R \
        ${per_chrom_mode_arg} \
        -s ${sample_id} \
        -t ${tumor_id} \
        -b ${baf_file} \
        -l ${logr_file} \
        ${position_scale_arg}\
        -o ./ \
        -r ${plot_resolution} \
        ${reference_dict_arg}
    """
}
