// Define the `FLAGSTAT` process that aligns stats
process FLAGSTAT {
    container = 'glebusasha/bwa_samtools'
    tag "$bamFile"
    cpus params.cpus
    publishDir "${params.outdir}/${workflow.start}[${workflow.runName}]/FLAGSTAT"
//	  debug true
//    errorStrategy 'ignore'

    input:
    path bamFile
    
    output:
    path "*.flagstat", emit: flagstat
    
    script:
    """
    samtools flagstat $bamFile > ${bamFile.baseName}.flagstat
    """
}
