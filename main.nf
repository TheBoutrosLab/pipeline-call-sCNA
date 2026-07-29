#!/usr/bin/env nextflow

nextflow.enable.dsl=2

include { run_validate_PipeVal } from './external/pipeline-Nextflow-module/modules/PipeVal/validate/main.nf'
include { indexFile } from './external/pipeline-Nextflow-module/modules/common/indexFile/main.nf'
include { workflow_cnv_facets } from './module/workflow-cnv_facets.nf'
include { workflow_battenberg } from './module/workflow-battenberg.nf'

log.info """\
=====================================
C A L L - S C N A  P I P E L I N E
=====================================
Boutros Lab

Current Configuration:

    - pipeline:
        name: ${workflow.manifest.name}
        version: ${workflow.manifest.version}

    - input:
        samples: ${params.samples_to_process}
        genome build: ${params.genome_build}
        target regions: ${params.target_bed}
        algorithms: ${params.algorithm}
        sample_sex: "${params.sample_sex}"
        battenberg_reference: "${params.battenberg_reference}"
        dbSNP_file: "${params.dbSNP_file}"
        reference_dict: "${params.reference_dict}"
        position_scale: "${params.position_scale}"

    - output:
        output_dir_base: ${params.output_dir_base}
        output: ${params.output_dir}
        log_output_dir: ${params.log_output_dir}

    - options:
        save_intermediate_files: ${params.save_intermediate_files}

    Tools Used:
        tool Battenberg: ${params.docker_image_battenberg}
        tool cnv_facets: ${params.docker_image_cnv_facets}
        tool BoutrosLabPlottingGeneral: ${params.docker_image_bpg}
        tool PipeVal: ${params.docker_image_validate}
        tool SigProfilerExtractor: ${params.docker_image_sigprofilerextractor}

    All parameters:
        ${params}

------------------------------------
Starting workflow...
------------------------------------
"""

workflow {
    /**
    *   Input channel processing
    */
    Channel.from(params.samples_to_process)
        .map{ sample -> ['index': indexFile(sample.path)] + sample }
        .set{ input_ch_samples_with_index }

    input_ch_samples_with_index
        .filter{ sample -> sample.sample_type == 'tumor' }
        .map{ sample -> [sample.id, sample.path, sample.index] }
        .set{ input_ch_tumor }

    input_ch_samples_with_index
        .filter{ sample -> sample.sample_type == 'normal' }
        .map{ sample -> [sample.id, sample.path, sample.index] }
        .set{ input_ch_normal }

    input_ch_tumor
        .combine(input_ch_normal)
        .set{ input_ch_paired_bams }

    Channel.fromPath(params.dbSNP_file, checkIfExists: true)
        .map{ dbSNP_file -> [dbSNP_file, indexFile(dbSNP_file.toString())] }
        .set{ input_ch_dbSNP_file }

    input_ch_samples_with_index
        .map{ sample -> [sample.path, sample.index] }
        .flatten()
        // .mix(input_ch_dbSNP_file.flatten())
        .set{ input_ch_validate }

    base_meta = Channel.value([
        'log_output_dir': "${params.log_output_dir}/process-log",
        'output_dir': params.output_dir_base
    ])

    module_meta = base_meta.map{ base_m ->
        base_m + [
            'docker_image': params.docker_image_validate
        ]
    }

    cnv_facets_meta = base_meta.map{ base_m ->
        base_m + [
            'workflow_output_dir': "${base_m.output_dir}/CNV_FACETS-${params.cnv_facets_version}"
        ]
    }

    battenberg_meta = base_meta.map{ base_m ->
        base_m + [
            'workflow_output_dir': "${base_m.output_dir}/Battenberg-${params.battenberg_version}"
        ]
    }

    /**
    *   Input validation
    */
    run_validate_PipeVal(
        module_meta.combine(input_ch_validate)
    )

    run_validate_PipeVal.out.validation_result
        .collectFile(
            name: 'input_validation.txt',
            storeDir: "${params.output_dir_base}/validation"
        )

    /**
    *   Call somatic copy-number variants with CNV_FACETS
    */
    if ('cnv_facets' in params.algorithm) {
        workflow_cnv_facets(
            cnv_facets_meta,
            input_ch_paired_bams,
            input_ch_dbSNP_file
        )
    }

    /**
    *   Call somatic copy-number variants with Battenberg
    */
    if ('Battenberg' in params.algorithm) {
        workflow_battenberg(
            battenberg_meta,
            input_ch_paired_bams
        )
    }
}
