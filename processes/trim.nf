// Define the `TRIM` process that performs quality trimming and filtering of reads
process TRIM{
    tag "${sid}"
    cpus params.cpus
    publishDir "${params.outdir}/TRIM"
//	  debug true
//    errorStrategy 'ignore'

    input:
    tuple val(sid), path(reads)

    output:
    tuple val(sid), path(fq_1_trimmed), path(fq_2_trimmed), emit: trimmed_reads
    path '*.html', emit: html, optional: true
    path '*.json', emit: json, optional: true

    script:
    fq_1_trimmed = sid + '_R1_P.fastq.gz'
    fq_2_trimmed = sid + '_R2_P.fastq.gz'
    """
    fastp -q 20 -l 20 --trim_poly_g --thread ${task.cpus} \
    --in1 ${reads[0]} \
    --in2 ${reads[1]}\
    --out1 $fq_1_trimmed \
    --out2 $fq_2_trimmed \
    --html ${sid}.fastp_stats.html \
    --json ${sid}.fastp_stats.json
    """
}
