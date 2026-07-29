nextflow.enable.dsl=2

include { call_SubclonalCopyNumber_Battenberg } from "./call_SubclonalCopyNumber_Battenberg"
include { suggest_refit_Battenberg } from "./suggest_refit_Battenberg"
include { workflow_plot_BAFLogR } from "./workflow-plot-BPG"
include { prepare_CNASignaturesInput_Battenberg } from "./prepare_CNASignaturesInput_Battenberg"
include { extract_CNASignatures_SigProfilerExtractor } from "./extract_CNASignatures_SigProfilerExtractor"
include { generate_checksum_PipeVal as generate_sha512_Battenberg } from "../external/pipeline-Nextflow-module/modules/PipeVal/generate-checksum/main.nf"

workflow workflow_battenberg {
    take:
    META
    paired_input_channel

    main:
    call_SubclonalCopyNumber_Battenberg(
        META,
        params.sample_sex,
        paired_input_channel
    )

    input_ch_script = Channel.fromPath("${projectDir}/script", checkIfExists: true)
    input_ch_reference_dict = Channel.fromPath(params.reference_dict, checkIfExists: true)

    suggest_refit_Battenberg(
        META,
        call_SubclonalCopyNumber_Battenberg.out.subclones_cna,
        call_SubclonalCopyNumber_Battenberg.out.default_refit_suggestion,
        params.sample,
        params.min_rho,
        params.min_psi,
        params.max_psi,
        input_ch_script
    )

    workflow_plot_BAFLogR(
        META,
        paired_input_channel,
        call_SubclonalCopyNumber_Battenberg.out.tumor_normal_baf_logr_files,
        input_ch_script,
        input_ch_reference_dict
    )

    prepare_CNASignaturesInput_Battenberg(
        META,
        params.sample,
        call_SubclonalCopyNumber_Battenberg.out.subclones_cna,
        input_ch_script
    )

    extract_CNASignatures_SigProfilerExtractor(
        META,
        params.sample,
        "BATTENBERG",
        prepare_CNASignaturesInput_Battenberg.out.battenberg_signature_input,
        input_ch_script
    )

    checksum_meta = META.map{ base_m ->
        base_m + [
            'output_dir': "${base_m.workflow_output_dir}/output",
            'docker_image': params.docker_image_validate
        ]
    }

    generate_sha512_Battenberg(
        checksum_meta.combine(call_SubclonalCopyNumber_Battenberg.out.subclones_cna)
    )

    emit:
    subclones_cna = call_SubclonalCopyNumber_Battenberg.out.subclones_cna
    default_refit_suggestion = call_SubclonalCopyNumber_Battenberg.out.default_refit_suggestion
    tumor_normal_baf_logr_files = call_SubclonalCopyNumber_Battenberg.out.tumor_normal_baf_logr_files
    custom_refit_suggestion = suggest_refit_Battenberg.out.custom_refit_suggestion
}
