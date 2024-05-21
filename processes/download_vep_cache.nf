// Define the `DOWNLOAD_VEP_CACHE` process that download VEP cache
process DOWNLOAD_VEP_CACHE {
    container = 'mgibio/vep_helper-cwl:vep_105.0_v1'
    tag "$cache_dir"
    publishDir "${params.outdir}/${workflow.start}[${workflow.runName}]/DOWNLOAD_VEP_CACHE"
    cpus params.cpus
//	  debug true
    errorStrategy 'ignore'
	
    input:
    path cache_dir

    script:
    """
    INSTALL.pl \
    -c $cache_dir \
    -a cf \
    -s homo_sapiens \
    -y GRCh38
    wget https://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh38/clinvar.vcf.gz -o $cache_dir/clinvar.vcf.gz
    """
}