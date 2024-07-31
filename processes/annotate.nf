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
    touch ${sid}.vep
    touch ${sid}.vep.html
    """
}

