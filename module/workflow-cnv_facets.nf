nextflow.enable.dsl=2

include { call_cnv_facets } from "./cnv_facets"
include { generate_checksum_PipeVal as generate_sha512_cnv_facets } from "../external/pipeline-Nextflow-module/modules/PipeVal/generate-checksum/main.nf"

workflow cnv_facets {
    take:
    META
    paired_input_channel
    dbSNP_file_channel

    main:
    target_bed_ch = Channel.fromPath(params.target_bed)

    call_cnv_facets(
        META,
        paired_input_channel,
        dbSNP_file_channel,
        target_bed_ch
    )

    checksum_meta = META.map{ base_m ->
        base_m + [
            'output_dir': "${base_m.workflow_output_dir}/output",
            'docker_image': params.docker_image_validate
        ]
    }

    generate_sha512_cnv_facets(
        checksum_meta.combine(
            call_cnv_facets.out.cnv_facets_vcf_and_index.flatten()
        )
    )
}
