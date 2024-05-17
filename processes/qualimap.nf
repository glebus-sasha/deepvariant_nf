// Define the `QUALIMAP` process that aligns stats
process QUALIMAP {
    tag "$bamFile"
    cpus params.cpus
    publishDir "${params.outdir}/QUALIMAP"
//	  debug true
//    errorStrategy 'ignore'

    input:
    path bamFile
    
    output:
    path "*"
    
    script:
    """
    qualimap bamqc -bam $bamFile
    """
}