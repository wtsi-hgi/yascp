
process GUZIP_VCF {
    tag "${samplename}"
    label 'process_medium'



    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "/software/hgi/containers/mercury_scrna_deconvolution_62bd56a-2021-12-15-4d1ec9312485.sif"
        //// container "/software/hgi/containers/mercury_scrna_deconvolution_latest.img"
    } else {
        container "mercury/scrna_deconvolution:62bd56a"
    }


    input: 
        tuple val(samplename), path(genotypes)
      

    output:
        tuple val(samplename), path("${samplename}.vcf"), emit: souporcell_vcf

    script:
      """
        bcftools view ${genotypes} -O v -o ${samplename}.vcf
      """
}
