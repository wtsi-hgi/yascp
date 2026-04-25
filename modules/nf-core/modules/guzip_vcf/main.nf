
process GUZIP_VCF {
    tag "${samplename}"
    label 'process_medium'



    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "${params.scrna_deconvolution}"
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
