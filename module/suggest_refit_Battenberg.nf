#!/usr/bin/env nextflow

log.info """\
---------------------------------------------------------------
        S U G G E S T  R E F I T - B A T T E N B E R G
---------------------------------------------------------------
Docker Images:
- docker_image_bpg: ${params.docker_image_bpg}
"""

include { generate_standard_filename; sanitize_string } from '../external/pipeline-Nextflow-module/modules/common/generate_standardized_filename/main.nf'

process suggest_refit_Battenberg {
    container params.docker_image_bpg

    tag "${sample_id}"

    publishDir "${META.workflow_output_dir}/output",
        pattern: "*.txt",
        mode: "copy",
        saveAs: {
            "${output_filename}_${sanitize_string(file(it).getName().replace("${sample_id}_", ""))}"
        }

    ext log_dir_suffix: { "-${sample_id}" }

    input:
        val(META)
        path(subclones_cna)
        path(default_refit_suggestion)
        val(sample_id)
        val(min_rho)
        val(min_psi)
        val(max_psi)
        path(script_dir)

    output:
        path "*custom_refit_suggestion.txt", optional: true, emit: custom_refit_suggestion

    script:
    output_filename = generate_standard_filename(
        "Battenberg-${params.battenberg_version}",
        params.dataset_id,
        sample_id,
        [:]
    )

    """
    set -euo pipefail

    Rscript ${script_dir}/suggest_refit_Battenberg.R \
        --subclones-file ${subclones_cna} \
        --default-refit-suggestions ${default_refit_suggestion} \
        --tumor-id ${sample_id} \
        --min-rho ${min_rho} \
        --min-psi ${min_psi} \
        --max-psi ${max_psi} \
        --platform-gamma 1 \
        --output-dir . \
        --util-script ${script_dir}/refit_util.R
    """
}
