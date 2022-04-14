
process SUBSET_GENOTYPE {
    tag "${samplename}.${sample_subset_file}"
    label 'process_medium'
    publishDir "${params.outdir}/subset_genotype/", mode: "${params.copy_mode}", pattern: "${samplename}.${sample_subset_file}.subset.vcf.gz"


    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "/software/hgi/containers/wtsihgi-nf_yascp_htstools-1.0.sif"
    } else {
        container "mercury/wtsihgi-nf_yascp_htstools-1.0"
    }

    input:
    tuple val(samplename), path(donor_vcf),path(donor_vcf_csi), val(sample_subset_file)


    output:
    tuple val(samplename), path("${samplename}.subset.vcf.gz"),path("${samplename}.subset.vcf.gz.csi"), emit: samplename_subsetvcf

    script:
    """
        echo ${sample_subset_file}
        #tabix -p vcf ${donor_vcf} || echo 'not typical VCF'
        bcftools view ${donor_vcf} -s ${sample_subset_file} -Oz -o ${samplename}.subset.vcf.gz
        bcftools index ${samplename}.subset.vcf.gz
        rm ${donor_vcf}.tbi || echo 'not typical VCF'
    """
}
