
process GUZIP_VCF {
    tag "${samplename}"
    label 'process_medium'



    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "/software/hgi/containers/scrna_deconvolution_latest.img"
    } else {
        container "quay.io/biocontainers/multiqc:1.10.1--py_0"
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
