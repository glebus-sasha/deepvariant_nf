// Define the `ANNOTATE` process that performs annotation
process ANNOTATE {
    container = 'ensemblorg/ensembl-vep:latest'
    containerOptions "-B ${params.vepcache}:/opt/vep/.vep"
    tag "$vcf"
    publishDir "${params.outdir}/${workflow.start}[${workflow.runName}]/ANNOTATE"
    cpus params.cpus
    debug true
    cache "lenient"
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
    --fork ${task.cpus} \
    --cache \
    --everything \
    --species homo_sapiens \
    --custom file=/opt/vep/.vep/clinvar.vcf.gz,short_name=ClinVar,format=vcf,type=exact,coords=0,fields=CLNSIG%CLNREVSTAT%CLNDN \
    --offline \
    --assembly GRCh38 \
    --dir ${params.vepcache}
    """
}
