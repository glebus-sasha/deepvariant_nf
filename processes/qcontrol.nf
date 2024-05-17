// Define the `QCONTROL` process that performs quality trimming and filtering of reads
process QCONTROL{
    tag "${sid}"
    cpus params.cpus
    publishDir "${params.outdir}/QCONTROL"
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
