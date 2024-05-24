// Define the `BAMINDEX` process that prepares the bam file indices
process BAMINDEX {
    container = 'glebusasha/bwa_samtools'
    tag "$bamFile"
    publishDir "${params.outdir}/${workflow.start.format('yyyy-MM-dd_HH-mm-ss')}_${workflow.runName}/BAMINDEX"
    cpus params.cpus
//	  debug true
//    errorStrategy 'ignore'
    input:
    tuple val(sid), path(bamFile)

    output:
    tuple val(sid), path('*.bai'), path(bamFile), emit: bai

    script:
    """
	samtools index ${bamFile}
    """
}