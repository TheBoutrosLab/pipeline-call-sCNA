nextflow.enable.dsl=2

include { plot_BAFLogR_BPG } from './bpg_baf_logr.nf'

def isSampleTypeFile(file, sampletype) {
    // match filenames with pattern sampletype and end with BAF.tab or LogR.tab
    return file.toString() =~ sampletype+"BAF.tab" || file.toString() =~ sampletype+"LogR.tab"
}

workflow workflow_plot_BAFLogR {
    take:
    META
    paired_input_channel
    baf_logr_data_ch
    script_dir_ch
    reference_dict_ch

    main:
    tumor_id_ch = paired_input_channel.map{
        tumor_id, tumor_bam, tumor_bam_index, normal_id, normal_bam, normal_bam_index ->
            tumor_id
    }

    normal_id_ch = paired_input_channel.map{
        tumor_id, tumor_bam, tumor_bam_index, normal_id, normal_bam, normal_bam_index ->
            normal_id
    }

    normal_baf_logr = baf_logr_data_ch.flatten()
        .filter {
            isSampleTypeFile(it, params.sample+"_normal")
        }
        .collect()
        .map { it -> it.sort() }

    normal_baf_logr_ch = normal_id_ch.combine(normal_baf_logr)

    tumor_baf_logr = baf_logr_data_ch.flatten()
        .filter {
            isSampleTypeFile(it, params.sample+"_mutant")
        }
        .collect()
        .map { it -> it.sort() }

    tumor_baf_logr_ch = tumor_id_ch.combine(tumor_baf_logr)

    normal_tumor_baf_logr = normal_baf_logr_ch.concat(tumor_baf_logr_ch)

    per_chrom_mode_flags = ["enabled", "disabled"]

    plot_BAFLogR_BPG(
        META,
        params.position_scale,
        params.sample,
        normal_tumor_baf_logr,
        per_chrom_mode_flags,
        script_dir_ch,
        params.bpg_plot_resolution,
        reference_dict_ch
    )
}
