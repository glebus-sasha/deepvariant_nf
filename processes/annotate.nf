// Define the `ANNOTATE` process that performs annotation
process ANNOTATE {
    container = 'mgibio/vep_helper-cwl:vep_105.0_v1'
    tag "$vcf"
    publishDir "${params.outdir}/${workflow.start}[${workflow.runName}]/ANNOTATE"
    cpus params.cpus
	  debug true
//    errorStrategy 'ignore'
	
    input:
    path vcf

    output:
    path '*.vep.vcf', emit: vep
    path '*.vep.html', emit: html

    script:
    """
    vep \
    -i $vcf \
    -o ${vcf.baseName}.vep.vcf \
    --stats_file ${vcf.baseName}.vep.html \
    --cache \
    --dir ${params.vepcache} \
    --vcf \ 
    --fork $task.cpus \
    --everything \
    --custom file=${params.vepcache}/clinvar.vcf.gz,short_name=ClinVar,format=vcf,type=exact,coords=0,fields=CLNSIG%CLNREVSTAT%CLNDN
    """
}
