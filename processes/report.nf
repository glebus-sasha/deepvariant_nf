// Define the `REPORT` process that performs report
process REPORT {
    container = 'staphb/multiqc:latest'
    tag "$flagstat"
    publishDir "${params.outdir}/${workflow.start.format('yyyy-MM-dd_HH-mm-ss')}_${workflow.runName}/REPORT"
    cpus params.cpus
//	  debug true
//    errorStrategy 'ignore'
	
    input:
    path flagstat
    path qualimap
    path vep

    output:
    path '*.html', emit: html

    script:
    """
    multiqc $flagstat $qualimap $vep
    """
}