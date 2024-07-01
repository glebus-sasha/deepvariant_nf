// Define the `ALIGN` process that aligns reads to the reference genome
process ALIGN {
    container = 'glebusasha/bwa_samtools'
    tag "$reference ${sid} $bedfile"
    cpus params.cpus
    publishDir "${params.outdir}/${workflow.start.format('yyyy-MM-dd_HH-mm-ss')}_${workflow.runName}/ALIGN"
//	  debug true
//    errorStrategy 'ignore'

    input:
    tuple val(sid), path(reads)
    path reference
    path idx
    path bedfile
    
    output:
    tuple val(sid), path("*.sorted.bam"), emit: bam
    
    script:
    def bed_option = bedfile.getBaseName() == 'dummy' ? "" : "-L ${bedfile}"    // If the base name of bedfile is 'dummy', set bed_option to an empty string
    """
        bwa mem \
            -t ${task.cpus} ${reference} ${reads[0]} ${reads[1]} | \
        samtools view -bh ${bed_option} | \
        samtools sort -o ${sid}.sorted.bam

    """
}