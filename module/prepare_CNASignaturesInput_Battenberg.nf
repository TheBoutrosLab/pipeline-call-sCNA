#!/usr/bin/env nextflow

log.info """\
---------------------------------------------------------------------------------
    P R E P A R E  C N A  S I G N A T U R E S  I N P U T  B A T T E N B E R G
---------------------------------------------------------------------------------
Docker Images:
- docker_image_sigprofilerextractor: ${params.docker_image_sigprofilerextractor}
"""

include { generate_standard_filename; sanitize_string } from '../external/pipeline-Nextflow-module/modules/common/generate_standardized_filename/main.nf'

/**
* Prepare Battenberg subclonal copy-number output for CNA signature extraction.
*/
process prepare_CNASignaturesInput_Battenberg {
    container params.docker_image_sigprofilerextractor

    tag "${tumor_id}"

    publishDir "${META.workflow_output_dir}/QC",
        pattern: "*seg_input.txt",
        mode: "copy",
        saveAs: {
            "${output_filename}_${sanitize_string(file(it).getName().replace("${tumor_id}_", ""))}"
        }

    ext log_dir_suffix: { "-${tumor_id}" }

    input:
        val(META)
        val(tumor_id)
        path(battenberg_subclones_file)
        path(script_dir)

    output:
        path "${tumor_id}*seg_input.txt", emit: battenberg_signature_input

    script:
    output_filename = generate_standard_filename(
        "Battenberg-${params.battenberg_version}",
        params.dataset_id,
        tumor_id,
        [:]
    )

    """
    set -euo pipefail

    Rscript ${script_dir}/prep_input_battenberg.R \
        --tumor-id ${tumor_id} \
        --subclones-file ${battenberg_subclones_file} \
        --output-dir ./
    """
}
