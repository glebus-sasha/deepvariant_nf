// Define the `VARCALL` process that performs variant calling
process VARCALL {
    tag "$reference $bamFile"
    publishDir "${params.outdir}/VARCALL"
    cpus params.cpus
//    debug true
    errorStrategy 'ignore'
	
    input:
    path reference
    path bamFile
    path bai
    path fai

    output:
    path "${bamFile.baseName}.vcf.gz", emit:vcf
    path "${bamFile.baseName}.g.vcf.gz", emit: gvcf
    path '*.html', emit: html
    
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