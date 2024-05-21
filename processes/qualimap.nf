// Define the `QUALIMAP` process that aligns stats
process QUALIMAP {
    container = 'pegi3s/qualimap:2.2.1_ubuntu19.01'
    tag "$bamFile"
    cpus params.cpus
    publishDir "${params.outdir}/${workflow.start}[${workflow.runName}]/QUALIMAP"
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