#!/usr/bin/env nextflow

log.info """\
------------------------------------
        C N V  F A C E T S
------------------------------------
Docker Images:
- docker_image_cnv_facets: ${params.docker_image_cnv_facets}
"""

include { generate_standard_filename } from '../external/pipeline-Nextflow-module/modules/common/generate_standardized_filename/main.nf'

/**
* Call somatic copy-number variants with CNV_FACETS.
*/
process call_cnv_facets {
    container params.docker_image_cnv_facets

    tag "${tumor_id}-${normal_id}"

    publishDir "${META.workflow_output_dir}/output",
        pattern: "${output_filename}.*",
        mode: "copy"

    ext log_dir_suffix: { "-${tumor_id}-${normal_id}" }

    input:
        val(META)
        tuple val(tumor_id), path(tumor_bam), path(tumor_bam_index), \
            val(normal_id), path(normal_bam), path(normal_bam_index)
        tuple path(dbSNP_file), path(dbSNP_file_index)
        path(target_bed)

    output:
        path "${output_filename}.vcf.gz*", emit: cnv_facets_vcf_and_index
        path "${output_filename}.csv.gz"
        path "${output_filename}.spider.pdf"
        path "${output_filename}.cnv.png"

    script:
    output_filename = generate_standard_filename(
        "CNV_FACETS-${params.cnv_facets_version}",
        params.dataset_id,
        tumor_id,
        [:]
    )

    target_file_arg = params.use_target_bed ? "--targets ${target_bed}" : ""
    no_cov_plot_arg = params.no_cov_plot ? "--no-cov-plot" : ""
    snp_count_orphans_arg = params.snp_count_orphans ? "--snp-count-orphans" : ""

    """
    set -euo pipefail

    cnv_facets.R \
        --snp-nprocs ${task.cpus} \
        --snp-vcf ${dbSNP_file} \
        --snp-mapq ${params.snp_mapq} \
        --snp-baq ${params.snp_baq} \
        --depth ${params.depth.sort().join(' ')} \
        --cval ${params.cval.sort().join(' ')} \
        --nbhd-snp ${params.nbhdsnp} \
        --rnd-seed ${params.rnd_seed} \
        --gbuild ${params.genome_build_dict[params.genome_build]} \
        ${snp_count_orphans_arg} \
        ${no_cov_plot_arg} \
        ${target_file_arg} \
        --snp-tumour ${tumor_bam} \
        --snp-normal ${normal_bam} \
        --out ${output_filename}
    """
}
