// Define the `FLAGSTAT` process that aligns stats
process FLAGSTAT {
    container = 'glebusasha/bwa_samtools'
    tag "$bamFile"
    cpus params.cpus
    publishDir "${params.outdir}/${workflow.start.format('yyyy-MM-dd_HH-mm-ss')}_${workflow.runName}/FLAGSTAT"
//	  debug true
//    errorStrategy 'ignore'

    input:
    tuple val(sid), path(bamFile)
    
    output:
    path "*.flagstat", emit: flagstat
    
    script:
    """
    samtools flagstat $bamFile > ${sid}.flagstat
    """
}
