#!/usr/bin/env nextflow

// Include processes
include { REFINDEX } from './processes/refindex.nf'
include { QCONTROL } from './processes/qcontrol.nf'
include { TRIM } from './processes/trim.nf'
include { ALIGN } from './processes/align.nf'
include { FLAGSTAT } from './processes/flagstat.nf'
include { QUALIMAP } from './processes/qualimap.nf'
include { FAINDEX } from './processes/faindex.nf'
include { BAMINDEX } from './processes/bamindex.nf'
include { VARCALL } from './processes/varcall.nf'
include { ANNOTATE } from './processes/annotate.nf'
include { REPORT } from './processes/report.nf'

// Logging pipeline information
log.info """\
    ==========================================
     D E E P V A R I A N T   P I P E L I N E
    ==========================================

    reference: ${params.reference}
    reads    : ${params.reads}
    outdir   : ${params.outdir}
    """
    .stripIndent(true)

// Define help
if ( params.help ) {
    help = """main.nf: This repository contains a Nextflow pipeline for analyzing 
            |Next-Generation Sequencing (NGS) data using octopus 
            |
            |Required arguments:
            |   --reference     Location of the reference file.
            |                   [default: ${params.reference}]
            |   --reads         Location of the input file file.
            |                   [default: ${params.reads}]
            |   --outdir        Location of the output file file.
            |                   [default: ${params.outdir}]
            |
            |Optional arguments:
            |   -profile        <docker/singularity>
            |
""".stripMargin()
    // Print the help with the stripped margin and exit
    println(help)
    exit(0)
}

// Define the input channel for FASTQ files, if provided
input_fastqs = params.reads ? Channel.fromFilePairs(params.reads, checkIfExists: true) : null

// Define the input channel for bwa index files, if provided
bwaidx = params.bwaidx ? Channel.fromPath(params.bwaidx, checkIfExists: true).collect() : null


// Define the workflow
workflow {
    // Make the pipeline reports directory if it needs
    if (params.reports) {
        def pipeline_report_dir = new File("${params.outdir}/pipeline_info")
        pipeline_report_dir.mkdirs()
    }

    QCONTROL(input_fastqs)
    TRIM(input_fastqs)
    if( !params.prebuild ) {
        REFINDEX(params.reference)
        ALIGN(TRIM.out.trimmed_reads, params.reference, REFINDEX.out)
    }
    else {
        ALIGN(TRIM.out.trimmed_reads, params.reference, bwaidx)
    }
    FLAGSTAT(ALIGN.out.bam)
    QUALIMAP(ALIGN.out.bam)
    FAINDEX(params.reference)
    BAMINDEX(ALIGN.out.bam)
    VARCALL(params.reference, ALIGN.out.bam, BAMINDEX.out.bai, FAINDEX.out.fai)
    ANNOTATE(VARCALL.out.vcf)
    REPORT(TRIM.out.json.collect(), QCONTROL.out.zip.collect(), FLAGSTAT.out.flagstat.collect(), QUALIMAP.out.collect(), ANNOTATE.out.html.collect())
}

// Log pipeline execution summary on completion
workflow.onComplete {
    log.info """\
        Pipeline execution summary
        ---------------------------
        Completed at: ${workflow.complete}
        Duration    : ${workflow.duration}
        Success     : ${workflow.success}
        workDir     : ${workflow.workDir}
        exit status : ${workflow.exitStatus}
        """
        .stripIndent()
        
    log.info ( workflow.success ? "\nDone" : "\nOops" )
}




