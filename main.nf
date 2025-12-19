#!/usr/bin/env nextflow
/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    taylordgenes/variant_calling
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Github : https://github.com/taylordgenes/variant_calling
----------------------------------------------------------------------------------------
*/

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT FUNCTIONS / MODULES / SUBWORKFLOWS / WORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { VARIANT_CALLING  } from './workflows/variant_calling'
include { PIPELINE_INITIALISATION } from './subworkflows/local/utils_nfcore_variant_calling_pipeline'
include { PIPELINE_COMPLETION     } from './subworkflows/local/utils_nfcore_variant_calling_pipeline'
include { getGenomeAttribute      } from './subworkflows/local/utils_nfcore_variant_calling_pipeline'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    GENOME PARAMETER VALUES
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

// TODO nf-core: Remove this line if you don't need a FASTA file
//   This is an example of how to use getGenomeAttribute() to fetch parameters
//   from igenomes.config using `--genome`
params.fasta = getGenomeAttribute('fasta')
params.fasta_fai = getGenomeAttribute('fasta_fai')
params.bwamem2 = getGenomeAttribute('bwamem2')

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    NAMED WORKFLOWS FOR PIPELINE
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//
// WORKFLOW: Run main analysis pipeline depending on type of input
//
workflow TAYLORDGENES_VARIANT_CALLING {

    take:
    samplesheet // channel: samplesheet read in from --input
    main:

    // Wrap iGenomes references into channels with meta maps
    ch_fasta       = Channel.value([ [id: params.genome], file(params.fasta) ])
    ch_fasta_fai   = Channel.value([ [id: params.genome], file(params.fasta_fai) ])
    ch_bwamem2     = Channel.value([ [id: params.genome], file(params.bwamem2) ])

    // mosdepth
    ch_mosdepth_bed = params.mosdepth_bed ? Channel.fromPath(params.mosdepth_bed, checkIfExists: true)
        .map { bed -> [[id: bed.baseName], bed] }
        : Channel.value([[id:'null'], []])

    // verifybamid2 files
    ch_svd = Channel.value([
        file(params.verifybamid_ud,  checkIfExists: true),
        file(params.verifybamid_mu,  checkIfExists: true),
        file(params.verifybamid_bed, checkIfExists: true)
    ])

    // somalier files
    ch_somalier_sites     = Channel.value([ [id: 'sites'], file(params.somalier_sites, checkIfExists: true) ])
    ch_somalier_ref_files = Channel.fromPath("${params.somalier_ref_dir}/*.somalier", checkIfExists: true).collect()
    ch_labelled_somalier_files = ch_somalier_ref_files
        .map { ref_files -> [ [id:'somalier_ref'], file(params.somalier_labels, checkIfExists: true), ref_files ] }

        // VEP
    // ch_cache  = Channel.value( file(params.vep_cache) )
    // ch_extra  = Channel.value( file(params.vep_extra_dir) ) // plugins/customs folder

    // // VCFs (bgzipped + tabixed earlier)
    // ch_vcfs = input_vcf_paths
    // .map { v -> [[id: v.baseName], file(v)] }
    // // if you want to pass a per-sample extras dir, otherwise use Channel.value(null)
    // .map { meta_vcf -> tuple(meta_vcf[0], meta_vcf[1], file(params.vep_custom_bundle)) }


    //
    // WORKFLOW: Run pipeline
    //
    VARIANT_CALLING (
        samplesheet,
        ch_fasta,
        ch_fasta_fai,
        ch_bwamem2,
        ch_mosdepth_bed,
        ch_svd,
        ch_somalier_sites,
        ch_labelled_somalier_files
    )
    emit:
    multiqc_report = VARIANT_CALLING.out.multiqc_report // channel: /path/to/multiqc_report.html
}
/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow {

    main:
    //
    // SUBWORKFLOW: Run initialisation tasks
    //
    PIPELINE_INITIALISATION (
        params.version,
        params.validate_params,
        params.monochrome_logs,
        args,
        params.outdir,
        params.input,
        params.help,
        params.help_full,
        params.show_hidden
    )

    //
    // WORKFLOW: Run main workflow
    //
    TAYLORDGENES_VARIANT_CALLING (
        PIPELINE_INITIALISATION.out.samplesheet
    )
    //
    // SUBWORKFLOW: Run completion tasks
    //
    PIPELINE_COMPLETION (
        params.email,
        params.email_on_fail,
        params.plaintext_email,
        params.outdir,
        params.monochrome_logs,
        params.hook_url,
        TAYLORDGENES_VARIANT_CALLING.out.multiqc_report
    )
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
