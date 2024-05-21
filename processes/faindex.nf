// Define the `FAINDEX` process that prepares the reference genome indices
process FAINDEX {
    container = 'glebusasha/bwa_samtools'
    tag "$reference"
    publishDir "${params.outdir}/${workflow.start}[${workflow.runName}]/FAINDEX"
    cpus params.cpus
//	  debug true
//    errorStrategy 'ignore'
    input:
    path reference

    output:
    path '*.fai', emit: fai

    script:
    """
    samtools faidx "$reference"
    """
}