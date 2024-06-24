// Define the `VARCALL` process that performs variant calling
process VARCALL {
    container = 'google/deepvariant:1.6.1'
    tag "$reference $bamFile"
    publishDir "${params.outdir}/${workflow.start.format('yyyy-MM-dd_HH-mm-ss')}_${workflow.runName}/VARCALL"
    cpus params.cpus
    cache "lenient" 
//    debug true
    errorStrategy 'ignore'
	
    input:
    path reference
    tuple val(sid), path(bai), path(bamFile)
    path fai

    output:
    val sid
    tuple val(sid), path("${sid}.vcf.gz"),      emit: vcf
    path "${sid}.g.vcf.gz",                     emit: gvcf
    path '*.html',                              emit: html
    
    script:
    """
    /opt/deepvariant/bin/run_deepvariant \
    --model_type=WES \
    --ref=$reference \
    --reads=$bamFile \
    --output_vcf=${sid}.vcf.gz \
    --output_gvcf=${sid}.g.vcf.gz \
    --num_shards=${task.cpus} 
    """
}