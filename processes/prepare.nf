// Define the `PREPARE` process that prepares the reference genome indices
process PREPARE {
    tag "$bamFile $reference"
    publishDir "${params.outdir}/PREPARE"
    cpus params.cpus
//	  debug true
//    errorStrategy 'ignore'
    input:
    path reference
    path bamFile

    output:
    path '*.sorted.bam.bai', emit: bai
    path '*.fai', emit: fai

    script:
    """
	samtools index ${bamFile.baseName}.bam
    samtools faidx "$reference"
    """
}