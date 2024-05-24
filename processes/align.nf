// Define the `ALIGN` process that aligns reads to the reference genome
process ALIGN {
    container = 'glebusasha/bwa_samtools'
    tag "$reference ${sid}"
    cpus params.cpus
    publishDir "${params.outdir}/${workflow.start.format('yyyy-MM-dd_HH-mm-ss')}_${workflow.runName}/ALIGN"
//	  debug true
//    errorStrategy 'ignore'

    input:
    tuple val(sid), path(reads1), path(reads2)
    path reference
    path idx
    
    output:
    tuple val(sid), path("*.sorted.bam"), emit: bam
    
    script:
    """
    bwa mem \
    -t ${task.cpus} ${reference} ${reads1} ${reads2} | \
    samtools view -bh | \
    samtools sort -o ${sid}.sorted.bam
    """
}