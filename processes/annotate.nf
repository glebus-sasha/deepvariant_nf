// Define the `ANNOTATE` process that performs annotation
process ANNOTATE {
    container = 'ensemblorg/ensembl-vep:latest'
//    containerOptions "-B ${params.vepcache}:/opt/vep/.vep"
    tag "$vcf"
    publishDir "${params.outdir}/${workflow.start.format('yyyy-MM-dd_HH-mm-ss')}_${workflow.runName}/ANNOTATE"
//    debug true
    cache "lenient"
    errorStrategy 'ignore'

    input:
    tuple val(sid), path(vcf)
    path vep_cache
    path reference
    path clinvar_gz
    path clinvar_tbi

    output:
    path "${sid}.vep", emit: vep
    path "${sid}.vep.html", emit: html

    script:
    """
    vep \
    -i $vcf \
    -o ${sid}.vep \
    --fork ${task.cpus} \
    --cache \
    --dir_cache ${vep_cache} \
    --everything \
    --species homo_sapiens \
    --custom file=${clinvar_gz},short_name=ClinVar,format=vcf,type=exact,coords=0,fields=CLNSIG%CLNREVSTAT%CLNDN \
    --offline \
    --assembly GRCh38
    touch ${sid}.vep.html
    """
}

