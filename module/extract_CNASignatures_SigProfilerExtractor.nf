#!/usr/bin/env nextflow

log.info """\
--------------------------------------
        C N A  S I G N A T U R E S
--------------------------------------
Docker Images:
- docker_image_sigprofilerextractor: ${params.docker_image_sigprofilerextractor}
"""

include { generate_standard_filename; sanitize_string } from '../external/pipeline-Nextflow-module/modules/common/generate_standardized_filename/main.nf'

def boolToFlag(paramValue, flagName) {
    return paramValue ? "--${flagName}" : ""
}

process extract_CNASignatures_SigProfilerExtractor {
    container params.docker_image_sigprofilerextractor

    tag "${tumor_id}-${cna_caller}"

    publishDir "${META.workflow_output_dir}/QC",
        mode: "copy",
        saveAs: {
            "${output_filename}_${sanitize_string(file(it).getName().replace("${tumor_id}_", ""))}"
        }

    ext log_dir_suffix: { "-${tumor_id}-${cna_caller}" }

    input:
        val(META)
        val(tumor_id)
        val(cna_caller)
        path(cna_signature_input)
        path(script_dir)

    output:
        path "BATTENBERG.CNV48.matrix.tsv", emit: cna_signature_matrix
        path "CNV48/", emit: cna_signatures
        path "CosmicTemplates/", emit: cosmic_templates
        path "JOB_METADATA.txt", emit: job_metadata
        path "Seeds.txt", emit: seeds

    script:
    output_filename = generate_standard_filename(
        "SigProfilerExtractor-${params.sigprofilerextractor_version}",
        params.dataset_id,
        tumor_id,
        [:]
    )

    """
    set -euo pipefail

    python ${script_dir}/run_sigprofilerextractor.py \
        --cna-caller ${cna_caller} \
        --sample ${tumor_id} \
        --input ${cna_signature_input} \
        --output-dir ./ \
        --reference-genome ${params.genome_build} \
        --opportunity-genome ${params.genome_build} \
        --context-type ${params.context_type} \
        ${boolToFlag(params.exome, "exome")} \
        --min-signatures ${params.min_max_signatures.sort()[0]} \
        --max-signatures ${params.min_max_signatures.sort()[1]} \
        --nmf-replicates ${params.nmf_replicates} \
        ${boolToFlag(params.resample, "resample")} \
        --seeds ${params.seeds} \
        --matrix-normalization ${params.matrix_normalization} \
        --nmf-init ${params.nmf_init} \
        --precision ${params.precision} \
        --min-nmf-iterations ${params.min_max_nmf_iterations.sort()[0]} \
        --max-nmf-iterations ${params.min_max_nmf_iterations.sort()[1]} \
        --nmf-test-conv ${params.nmf_test_cov} \
        --nmf-tolerance ${params.nmf_tolerance} \
        --cpu ${task.cpus} \
        ${boolToFlag(params.gpu, "gpu")} \
        --batch-size ${params.batch_size} \
        --stability ${params.stability} \
        --min-stability ${params.min_stability} \
        --combined-stability ${params.combined_stability} \
        ${boolToFlag(params.allow_stability_drop, "allow-stability-drop")} \
        --cosmic-version ${params.cosmic_version} \
        ${boolToFlag(params.make_decomposition_plots, "make-decomposition-plots")} \
        ${boolToFlag(params.collapse_to_SBS96, "collapse-to-SBS96")} \
        ${boolToFlag(params.get_all_signature_matrices, "get-all-signature-matrices")} \
        ${boolToFlag(params.export_probabilities, "export-probabilities")} \
        --volume ./
    """
}
