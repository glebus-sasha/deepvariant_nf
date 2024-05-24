// Define the `QUALIMAP` process that aligns stats
process QUALIMAP {
    container = 'pegi3s/qualimap:2.2.1_ubuntu19.01'
    tag "$bamFile"
    cpus params.cpus
    publishDir "${params.outdir}/${workflow.start.format('yyyy-MM-dd_HH-mm-ss')}_${workflow.runName}/QUALIMAP"
//	  debug true
//    errorStrategy 'ignore'

    input:
    tuple val(sid), path(bamFile)
    
    output:
    path "*"
    
    script:
    """
    qualimap bamqc -bam $bamFile
    """
}