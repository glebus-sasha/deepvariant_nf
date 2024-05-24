// Define the `VARCALL` process that performs variant calling
process VARCALL {
    container = 'google/deepvariant:1.6.1'
    tag "$reference $bamFile"
    publishDir "${params.outdir}/${workflow.start}[${workflow.runName}]/VARCALL"
    cpus params.cpus
    cache "lenient" 
//    debug true
    errorStrategy 'ignore'
	
    input:
    path reference
    path bamFile
    path bai
    path fai

    output:
    path "${bamFile.baseName}.vcf.gz",      emit:vcf
    path "${bamFile.baseName}.g.vcf.gz",    emit: gvcf
    path '*.html',                          emit: html
    
    script:
    """
    /opt/deepvariant/bin/run_deepvariant \
    --model_type=WES \
    --ref=$reference \
    --reads=$bamFile \
    --output_vcf=${bamFile.baseName}.vcf.gz \
    --output_gvcf=${bamFile.baseName}.g.vcf.gz \
    --num_shards=${task.cpus} 
    """
}