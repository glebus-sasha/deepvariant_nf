// Define the `ALIGN` process that aligns reads to the reference genome
process ALIGN {
    tag "$reference ${sid}"
    cpus params.cpus
    publishDir "${params.outdir}/ALIGN"
//	  debug true
//    errorStrategy 'ignore'

    input:
    tuple val(sid), path(reads1), path(reads2)
    path reference
    path idx
    
    output:
    path "${sid}.sorted.bam", emit: bam
    
    script:
    """
    bwa mem \
    -t ${task.cpus} ${reference} ${reads1} ${reads2} | \
    samtools view -bh | \
    samtools sort -o ${sid}.sorted.bam
    """
}