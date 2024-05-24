// Define the `QCONTROL` process that performs quality trimming and filtering of reads
process QCONTROL{
    container = 'staphb/fastqc:0.12.1'
    tag "${sid}"
    cpus params.cpus
    publishDir "${params.outdir}/${workflow.start.format('yyyy-MM-dd_HH-mm-ss')}_${workflow.runName}/QCONTROL", pattern: '*.html'
//	  debug true
//    errorStrategy 'ignore'

    input:
    tuple val(sid), path(reads)

    output:
    path "*.html", emit: html
    path "*.zip", emit: zip

    """
    fastqc $reads
    """
}
