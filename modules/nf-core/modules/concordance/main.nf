process CONCORDANCE_CALCLULATIONS {

    tag "${pool_id}"
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
        tuple val(pool_id), path("cell_concordance_table.tsv"), emit: concordances

    script:

        """
            echo ${pool_id}
            bcftools view ${vcf_exp} -R ${cell_vcf} -Oz -o sub_${pool_id}_Expected.vcf.gz
            if [ "${vcf_gt_match}" != "fake_file.fq" ]; then
                bcftools view ${vcf_gt_match} -R ${cell_vcf} -Oz -o sub_${pool_id}_GT_Matched.vcf.gz
            fi
            
            bcftools view  -i 'FORMAT/DP > 3' ${cell_vcf} -Oz -o sub_${pool_id}_cellSNP_dp_filter.vcf.gz
            concordance_calculations_donor_exclusive_dp.py --cpus $task.cpus --cell_vcf ${cell_vcf} --cell_vcf_dp sub_${pool_id}_cellSNP_dp_filter.vcf.gz --donor_assignments ${donor_table} --gt_match_vcf sub_${pool_id}_GT_Matched.vcf.gz --expected_vcf sub_${pool_id}_Expected.vcf.gz --cell_assignments ${cell_assignments}
        """
}


process COMBINE_FILES{

    tag "${pool_id}"
    label 'process_low'
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
        path(subsampling_table), 
        path(concordance_table))

    output:
        path("*.png")
        tuple val(pool_id), path("*joined_df_for_plots.tsv"), emit: joined_df_for_plots
        path("*joined_df_for_plots.tsv"), emit: file_joined_df_for_plots

    script:

        """
           combine_concordance.py -cc ${concordance_table} -sq ${subsampling_table} -name ${pool_id}
        """

}


process PLOT_CONCORDANCES_ALL{

    tag "${pool_id}"
    label 'process_low'
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "https://yascp.cog.sanger.ac.uk/public/singularity_images/mercury_scrna_deconvolution_62bd56a-2021-12-15-4d1ec9312485.sif"
        //// container "https://yascp.cog.sanger.ac.uk/public/singularity_images/mercury_scrna_deconvolution_latest.img"
    } else {
        container "mercury/scrna_deconvolution:62bd56a"
    }

    publishDir  path: "${params.outdir}/concordances",
                mode: "${params.copy_mode}",
                overwrite: "true"

    output:
        path("*.png")

    input:
        path(input_file_all)

    script:

        """
            combined_concordance_plots.py -cc ${input_file_all} -name all
        """


}
