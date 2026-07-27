#!/usr/bin/env nextflow

log.info """\
------------------------------------
        B A T T E N B E R G
------------------------------------
Docker Images:
- docker_image_battenberg: ${params.docker_image_battenberg}
"""

include { generate_standard_filename; sanitize_string } from '../external/pipeline-Nextflow-module/modules/common/generate_standardized_filename/main.nf'

/**
* Call somatic subclonal copy-number variants with Battenberg.
*/
process call_SubclonalCopyNumber_Battenberg {
    container params.docker_image_battenberg
    containerOptions "${params.container_mount_flag} ${params.battenberg_reference}:/opt/battenberg_reference/"

    tag "${tumor_id}-${normal_id}"

    publishDir "${META.workflow_output_dir}/output",
        pattern: "*.{txt,tab,png}",
        mode: "copy",
        saveAs: {
            !(file(it).getName().endsWith("${output_filename}_subclones.txt")) ? \
                "${output_filename}_${sanitize_string(file(it).getName().replace("${tumor_id}_", ""))}" : \
                    "${output_filename}_subclones.txt"
        }

    ext log_dir_suffix: { "-${tumor_id}-${normal_id}" }

    input:
        val(META)
        val(sample_sex)
        tuple val(tumor_id), path(tumor_bam), path(tumor_bam_index), \
            val(normal_id), path(normal_bam), path(normal_bam_index)

    output:
        path "${output_filename}_subclones.txt", emit: subclones_cna
        path "${tumor_id}_refit_suggestion.txt", emit: default_refit_suggestion
        path "${tumor_id}*.txt"
        path "${tumor_id}*.png"
        path "${tumor_id}*.tab", emit: tumor_normal_baf_logr_files
        path "sample_g.txt"

    script:
    output_filename = generate_standard_filename(
        "Battenberg-${params.battenberg_version}",
        params.dataset_id,
        tumor_id,
        [:]
    )

    """
    set -euo pipefail

    Rscript /usr/local/lib/R/site-library/Battenberg/example/battenberg_wgs.R \
        -t ${tumor_id} \
        -n ${normal_id} \
        --tb `readlink -f ${tumor_bam}` \
        --nb `readlink -f ${normal_bam}` \
        -o ./ \
        --sex ${sample_sex} \
        --cpu ${task.cpus} \
        --min_rho ${params.min_rho} \
        --max_rho ${params.max_rho} \
        --min_ploidy ${params.min_psi} \
        --max_ploidy ${params.max_psi} \
        --min_goodness_of_fit ${params.min_goodness_of_fit} \
        --balanced_threshold ${params.balanced_threshold} \
        --min_normal_depth ${params.min_normal_depth} \
        --min_base_qual ${params.min_base_qual} \
        --min_map_qual ${params.min_map_qual}

    mv ./${tumor_id}_subclones.txt ./${output_filename}_subclones.txt
    """
}
