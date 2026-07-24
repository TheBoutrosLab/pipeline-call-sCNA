#!/usr/bin/env nextflow

nextflow.enable.dsl=2

include { run_validate_PipeVal } from './external/pipeline-Nextflow-module/modules/PipeVal/validate/main.nf'
include { indexFile } from './external/pipeline-Nextflow-module/modules/common/indexFile/main.nf'

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
        dnSNP_file: "${params.dbSNP_file}"
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
        .map{ sample -> [sample.path, sample.index] }
        .flatten()
        .set{ input_ch_validate }

    base_meta = Channel.value([
        'log_output_dir': params.log_output_dir,
        'output_dir': params.output_dir_base
    ])

    module_meta = base_meta.map{ base_m ->
        base_m + [
            'log_output_dir': "${base_m.log_output_dir}/process-log"
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
}
