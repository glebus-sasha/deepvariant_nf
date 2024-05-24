// Define the `REFINDEX` process that creates the index of the genome
process REFINDEX {
    container = 'glebusasha/bwa_samtools'
    tag "$reference"
    publishDir "${params.outdir}/${workflow.start.format('yyyy-MM-dd_HH-mm-ss')}_${workflow.runName}/REFINDEX"
    cpus params.cpus
//	  debug true
//    errorStrategy 'ignore'

    input:
    path reference
    debug true

    output:
    path "*", emit: bwaidx

    script:
    """
    bwa index $reference
    """
}
