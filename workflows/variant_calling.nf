/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
include { FASTQC                 } from '../modules/nf-core/fastqc/main'
include { MULTIQC                } from '../modules/nf-core/multiqc/main'
include { FASTP                  } from '../modules/nf-core/fastp/main'
include { BWAMEM2_MEM            } from '../modules/nf-core/bwamem2/mem/main'
include { PICARD_COLLECTWGSMETRICS} from '../modules/nf-core/picard/collectwgsmetrics/main'
include { PICARD_MARKDUPLICATES  } from '../modules/nf-core/picard/markduplicates/main'
include { paramsSummaryMap       } from 'plugin/nf-schema'
include { paramsSummaryMultiqc   } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { softwareVersionsToYAML } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { methodsDescriptionText } from '../subworkflows/local/utils_nfcore_variant_calling_pipeline'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow VARIANT_CALLING {

    take:
    ch_samplesheet // channel: samplesheet read in from --input
    ch_fasta
    ch_fasta_fai
    ch_bwamem2
    main:

    ch_versions = Channel.empty()
    ch_multiqc_files = Channel.empty()

    //
    // MODULE: Run FastQC
    //
    FASTQC (
        ch_samplesheet
    )
    ch_multiqc_files = ch_multiqc_files.mix(FASTQC.out.zip.collect{it[1]})
    ch_versions = ch_versions.mix(FASTQC.out.versions.first())

    //
    // MODULE: Run Fastp
    //
    FASTP (
        ch_samplesheet,
        params.discard_trimmed_pass,
        params.save_trimmed_fail,
        params.save_merged
    )
    ch_multiqc_files = ch_multiqc_files.mix(FASTP.out.html.collect{ it[1] })
    ch_multiqc_files = ch_multiqc_files.mix(FASTP.out.json.collect{ it[1] })
    ch_versions = ch_versions.mix(FASTP.out.versions.first())

    //
    // MODULE: Run BWA-MEM2
    //
    BWAMEM2_MEM (
        FASTP.out.reads,
        ch_bwamem2,
        ch_fasta,
        params.sort_bam
    )
    ch_versions = ch_versions.mix(FASTP.out.versions.first())

    //
    // MODULE: Run Picard
    //
    PICARD_MARKDUPLICATES (
        BWAMEM2_MEM.out.bam,
        ch_fasta,
        ch_fasta_fai
    )
    ch_multiqc_files = ch_multiqc_files.mix(PICARD_MARKDUPLICATES.out.metrics.collect{ it[1] })

    intervallist = params.intervallist ? Channel.fromPath(params.intervallist).collect() : Channel.empty()
    ch_bam = BWAMEM2_MEM.out.bam            // tuple: (meta, bam)
    ch_bai = BWAMEM2_MEM.out.bai            // tuple: (meta, bai)
    ch_bam_bai = ch_bam.join(ch_bai, by: 0) // tuple: (meta, bam, bai)
    PICARD_COLLECTWGSMETRICS (
        ch_bam_bai,
        ch_fasta,
        ch_fasta_fai,
        intervallist
    )
    ch_multiqc_files = ch_multiqc_files.mix(PICARD_COLLECTWGSMETRICS.out.metrics.collect{ it[1] })
    ch_versions = ch_versions.mix(FASTP.out.versions.first())

    //
    // Collate and save software versions
    //
    softwareVersionsToYAML(ch_versions)
        .collectFile(
            storeDir: "${params.outdir}/pipeline_info",
            name:  'variant_calling_software_'  + 'mqc_'  + 'versions.yml',
            sort: true,
            newLine: true
        ).set { ch_collated_versions }


    //
    // MODULE: MultiQC
    //
    ch_multiqc_config        = Channel.fromPath(
        "$projectDir/assets/multiqc_config.yml", checkIfExists: true)
    ch_multiqc_custom_config = params.multiqc_config ?
        Channel.fromPath(params.multiqc_config, checkIfExists: true) :
        Channel.empty()
    ch_multiqc_logo          = params.multiqc_logo ?
        Channel.fromPath(params.multiqc_logo, checkIfExists: true) :
        Channel.empty()

    summary_params      = paramsSummaryMap(
        workflow, parameters_schema: "nextflow_schema.json")
    ch_workflow_summary = Channel.value(paramsSummaryMultiqc(summary_params))
    ch_multiqc_files = ch_multiqc_files.mix(
        ch_workflow_summary.collectFile(name: 'workflow_summary_mqc.yaml'))
    ch_multiqc_custom_methods_description = params.multiqc_methods_description ?
        file(params.multiqc_methods_description, checkIfExists: true) :
        file("$projectDir/assets/methods_description_template.yml", checkIfExists: true)
    ch_methods_description                = Channel.value(
        methodsDescriptionText(ch_multiqc_custom_methods_description))

    ch_multiqc_files = ch_multiqc_files.mix(ch_collated_versions)
    ch_multiqc_files = ch_multiqc_files.mix(
        ch_methods_description.collectFile(
            name: 'methods_description_mqc.yaml',
            sort: true
        )
    )

    MULTIQC (
        ch_multiqc_files.collect(),
        ch_multiqc_config.toList(),
        ch_multiqc_custom_config.toList(),
        ch_multiqc_logo.toList(),
        [],
        []
    )

    emit:multiqc_report = MULTIQC.out.report.toList() // channel: /path/to/multiqc_report.html
    versions       = ch_versions                 // channel: [ path(versions.yml) ]

}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
