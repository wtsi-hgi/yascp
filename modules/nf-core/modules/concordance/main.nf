process CONCORDANCE_CALCLULATIONS {

    tag "${samplename}"
    label 'process_medium'
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "https://yascp.cog.sanger.ac.uk/public/singularity_images/mercury_scrna_deconvolution_62bd56a-2021-12-15-4d1ec9312485.sif"
        //// container "https://yascp.cog.sanger.ac.uk/public/singularity_images/mercury_scrna_deconvolution_latest.img"
    } else {
        container "mercury/scrna_deconvolution:62bd56a"
    }

    publishDir  path: "${params.outdir}/concordances/${pool_id}",
                mode: "${params.copy_mode}",
                overwrite: "true"

    input:
        tuple(val(pool_id), 
        path(vcf_gt_match), 
        path(vcf_gt_match_csi),
        path(vcf_exp), 
        path(vcf_exp_csi),
        path(cell_vcf),
        path(donor_table),path(cell_assignments))

    output:
        path("cell_concordance_table.tsv", emit: concordances)

    script:

        """
            echo ${pool_id}
            bcftools view ${vcf_exp} -R ${cell_vcf} -Oz -o sub_${pool_id}_Expected.vcf.gz
            bcftools view ${vcf_gt_match} -R ${cell_vcf} -Oz -o sub_${pool_id}_GT_Matched.vcf.gz
            bcftools view  -i 'FORMAT/DP > 3' ${cell_vcf} -Oz -o sub_${pool_id}_cellSNP_dp_filter.vcf.gz
            concordance_calculations_donor_exclusive_dp.py --cpus $task.cpus --cell_vcf ${cell_vcf} --cell_vcf_dp sub_${pool_id}_cellSNP_dp_filter.vcf.gz --donor_assignments ${donor_table} --gt_match_vcf sub_${pool_id}_GT_Matched.vcf.gz --expected_vcf sub_${pool_id}_Expected.vcf.gz --cell_assignments ${cell_assignments}
        """
}
