// Define the `REPORT` process that performs report
process REPORT {
    container = 'staphb/multiqc:latest'
    tag "$flagstat"
    publishDir "${params.outdir}/${workflow.start}[${workflow.runName}]/REPORT"
    cpus params.cpus
//	  debug true
//    errorStrategy 'ignore'
	
    input:
    path fastp
    path fastqc
    path flagstat
    path qualimap
    path vep

    output:
    path '*.html', emit: html

    script:
    """
    multiqc $fastqc $fastp $flagstat $qualimap $vep
    """
}