// Define the `BAMINDEX` process that prepares the bam file indices
process BAMINDEX {
    container = 'glebusasha/bwa_samtools'
    tag "$bamFile"
    publishDir "${params.outdir}/${workflow.start}[${workflow.runName}]/BAMINDEX"
    cpus params.cpus
//	  debug true
//    errorStrategy 'ignore'
    input:
    path bamFile

    output:
    path '*.sorted.bam.bai', emit: bai

    script:
    """
	samtools index ${bamFile.baseName}.bam
    """
}