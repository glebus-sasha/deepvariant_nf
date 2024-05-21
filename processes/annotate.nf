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
    --custom ${params.vepcache}/clinvar.vcf.gz,ClinVar,vcf,exact,0,CLNSIG,CLNREVSTAT,CLNDN \
    --everything 
    """
}
