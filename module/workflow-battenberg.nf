nextflow.enable.dsl=2

include { call_SubclonalCopyNumber_Battenberg } from "./call_SubclonalCopyNumber_Battenberg"
include { generate_checksum_PipeVal as generate_sha512_Battenberg } from "../external/pipeline-Nextflow-module/modules/PipeVal/generate-checksum/main.nf"

workflow workflow_battenberg {
    take:
    META
    paired_input_channel

    main:
    sample_sex_ch = Channel.value(params.sample_sex)

    call_SubclonalCopyNumber_Battenberg(
        META,
        sample_sex_ch,
        paired_input_channel
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
}
