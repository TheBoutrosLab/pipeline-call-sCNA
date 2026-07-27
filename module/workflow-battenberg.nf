nextflow.enable.dsl=2

workflow workflow_battenberg {
    take:
    META
    paired_input_channel

    main:
    /*
    * Bundle the workflow metadata with each tumor-normal pair. Battenberg
    * processes will consume this channel as they are added to the workflow.
    */
    battenberg_input_channel = META.combine(paired_input_channel)

    emit:
    battenberg_input_channel
}
