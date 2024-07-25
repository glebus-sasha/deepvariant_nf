// Define the `VARCALL` process that performs variant calling
process VARCALL {
    container = 'google/deepvariant:1.6.1'
    tag "$reference $bamFile"
    publishDir "${params.outdir}/${workflow.start.format('yyyy-MM-dd_HH-mm-ss')}_${workflow.runName}/VARCALL"
//    cpus 1
//    memory '2 GB'
//    cache "lenient" 
//    debug true
    errorStrategy 'ignore'
	
    input:
    path reference
    tuple val(sid), path(bai), path(bamFile)
    path fai
    path bedfile

    output:
    val sid
    tuple val(sid), path("${sid}.vcf.gz"),      emit: vcf
    path "${sid}.g.vcf.gz",                     emit: gvcf
    path '*.html',                              emit: html
    
    script:
    def bed_option = bedfile.getBaseName() == 'dummy' ? "" : "--regions ${bedfile}"    // If the base name of bedfile is 'dummy', set bed_option to an empty string
    """
    /opt/deepvariant/bin/run_deepvariant \
    --model_type=WES \
    --ref=$reference \
    --reads=$bamFile \
    --output_vcf=${sid}.vcf.gz \
    --output_gvcf=${sid}.g.vcf.gz ${bed_option}\
    --num_shards=1
    """
}