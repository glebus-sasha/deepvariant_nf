#!/usr/bin/env nextflow

log.info """\
    ==========================================
     D E E P   V A R I A N T   P I P E L I N E
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

// Define the `REFINDEX` process that creates the index of the genome
process REFINDEX {
    tag "$reference"
//    publishDir "${params.outdir}/REFINDEX"

    input:
    path reference
    debug true

    output:
    path "*"

    script:
    """
    bwa index $reference
    """
}

// Define the `QCONTROL` process that performs quality trimming and filtering of reads
process QCONTROL{
    tag "${sid}"
    cpus params.cpus
    publishDir "${params.outdir}/QCONTROL"

    input:
    tuple val(sid), path(reads)

    output:
    tuple val(sid), path(fq_1_trimmed), path(fq_2_trimmed), emit: trimmed_reads
            file("${sid}.fastp_stats.html")

    script:
    fq_1_trimmed = sid + '_R1_P.fastq.gz'
    fq_2_trimmed = sid + '_R2_P.fastq.gz'
    """
    fastp -q 20 -l 140 --trim_poly_g --thread ${task.cpus} \
    --in1 ${reads[0]} \
    --in2 ${reads[1]}\
    --out1 $fq_1_trimmed \
    --out2 $fq_2_trimmed \
    --html ${sid}.fastp_stats.html
    """
}

// Define the `ALIGN` process that aligns reads to the reference genome
process ALIGN {
    tag "$reference ${sid}"
    cpus params.cpus
//    publishDir "${params.outdir}/ALIGN"
    debug true

    input:
    tuple val(sid), path(reads1), path(reads2)
    path reference
    path idx
    
    output:
    path "${sid}.sorted.bam"
    script:
    """
    bwa mem \
    -R "@RG\\tID:S41\\tSM:H1_U5\\tLB:M4\\tPU:Illumina" \
    -t ${task.cpus} ${reference} ${reads1} ${reads2} | \
    samtools view -bh | \
    samtools sort -o ${sid}.sorted.bam
    """
}

// Define the `PREPARE` process that prepares the reference genome indices
process PREPARE {
    tag "$bamFile $reference"
//    publishDir "${params.outdir}/PREPARE"
	
    input:
    path reference
    path bamFile

    output:
    file '*.sorted.bam.bai'
    file '*.fai'

    script:
    """
	samtools index ${bamFile.baseName}.bam
    samtools faidx $reference
    """
}

// Define the `VARCALL` process that performs variant calling
process VARCALL {
    tag "$reference $bamFile"
    publishDir "${params.outdir}/VARCALL"
    cpus params.cpus
	
    input:
    path reference
    path bamFile
    path bai
    path fai

    output:
    file "${bamFile.baseName}.vcf.gz"
    file "${bamFile.baseName}.g.vcf.gz"
    file '*.html'
    
    script:
    """
    /opt/deepvariant/bin/run_deepvariant \
    --model_type=WGS \
    --ref=$reference \
    --reads=$bamFile \
    --output_vcf=${bamFile.baseName}.vcf.gz \
    --output_gvcf=${bamFile.baseName}.g.vcf.gz \
    --num_shards=${task.cpus}
    """
}

// Define the `ANNOTATE` process that performs annotation
process ANNOTATE {
    tag "$vcf"
    publishDir "${params.outdir}/ANNOTATE"
	debug true
	
    input:
    path vcf

    output:
    file '*.vep.vcf'

    script:
    """
    vep --database -i $vcf -o ${vcf.baseName}.vep.vcf --everything
    """
}

// Define the input channels for FASTQ files, if provided
input_fastqs = params.reads ? Channel.fromFilePairs(params.reads, checkIfExists: true) : null

// Define the workflow
workflow {
		REFINDEX(params.reference)
		QCONTROL(input_fastqs)
       	ALIGN(QCONTROL.out[0], params.reference, REFINDEX.out)
      	PREPARE(params.reference, ALIGN.out)
       	VARCALL(params.reference, ALIGN.out, PREPARE.out[0], PREPARE.out[1])
        ANNOTATE(VARCALL.out[0])
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




