// Define the `ANNOTATE` process that performs annotation
process ANNOTATE {
    container = 'ensemblorg/ensembl-vep:latest'
    tag "$vcf"
    publishDir "${params.outdir}/${workflow.start}[${workflow.runName}]/ANNOTATE"
    cpus params.cpus
	  debug true
//    errorStrategy 'ignore'
	
    input:
    path vcf
    path gzi

    output:
    path '*.vep.vcf', emit: vep
    path '*.vep.html', emit: html

    script:
    """
    vep \
    -i $vcf \
    -o ${vcf.baseName}.vep.vcf \
    --stats_file ${vcf.baseName}.vep.html \
    --fork ${task.cpus} \
    --cache \
    --dir ${params.vepcache} \ 
    --custom file=${params.vepcache}/clinvar.vcf.gz,short_name=ClinVar,format=vcf,type=exact,coords=0,fields=CLNSIG%CLNREVSTAT%CLNDN
    --everything 
    """
}
