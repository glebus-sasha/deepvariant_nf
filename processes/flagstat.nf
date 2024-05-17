// Define the `FLAGSTAT` process that aligns stats
process FLAGSTAT {
    tag "$bamFile"
    cpus params.cpus
    publishDir "${params.outdir}/FLAGSTAT"
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
