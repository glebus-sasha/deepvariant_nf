// Define the `REFINDEX` process that creates the index of the genome
process REFINDEX {
    tag "$reference"
    publishDir "${params.outdir}/REFINDEX"
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
