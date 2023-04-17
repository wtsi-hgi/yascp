process CONCORDANCE_CALCLULATIONS {

    tag "${samplename}"
    label 'process_medium'
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "https://yascp.cog.sanger.ac.uk/public/singularity_images/wtsihgi_nf_scrna_qc_6bb6af5-2021-12-23-3270149cf265.sif"
    } else {
        container "wtsihgi/nf_scrna_qc:6bb6af5"
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
    // pool12, 
    // /lustre/scratch125/humgen/teams/hgi/mo11/oneK1k/extra0/work/23/b826d0054c4bb62450dc39d04a3f24/Genotype___GTMatchedSubset_pool12/Study_Merge_GTMatchedSubset_EUY8DDDZD_out.vcf.gz,
    //  /lustre/scratch125/humgen/teams/hgi/mo11/oneK1k/extra0/work/23/b826d0054c4bb62450dc39d04a3f24/Genotype___GTMatchedSubset_pool12/Study_Merge_GTMatchedSubset_EUY8DDDZD_out.vcf.gz.csi,
    //  /lustre/scratch125/humgen/teams/hgi/mo11/oneK1k/extra0/work/d2/34d2430908d3d27632fa0f8fb74357/Genotype___AllExpectedGT_pool12/Study_Merge_AllExpectedGT_SYIDTL7VN_out.vcf.gz,
    //  /lustre/scratch125/humgen/teams/hgi/mo11/oneK1k/extra0/work/d2/34d2430908d3d27632fa0f8fb74357/Genotype___AllExpectedGT_pool12/Study_Merge_AllExpectedGT_SYIDTL7VN_out.vcf.gz.csi, 
    // /lustre/scratch125/humgen/teams/hgi/mo11/oneK1k/extra0/work/da/3f5fd59e2498b18a3f9d64f256e187/stats_pool12_gt_donor_assignments.csv

    // pool53, 
    // /lustre/scratch125/humgen/teams/hgi/mo11/oneK1k/extra0/work/2b/f55757bbaa9287df4fa08a95c9572e/Genotype___GTMatchedSubset_pool53/Study_Merge_GTMatchedSubset_B9Y4IL5QE_out.vcf.gz, 
    // /lustre/scratch125/humgen/teams/hgi/mo11/oneK1k/extra0/work/2b/f55757bbaa9287df4fa08a95c9572e/Genotype___GTMatchedSubset_pool53/Study_Merge_GTMatchedSubset_B9Y4IL5QE_out.vcf.gz.csi,
    //  /lustre/scratch125/humgen/teams/hgi/mo11/oneK1k/extra0/work/ea/b902eb8a42b41fbeff7b25bd217bc3/Genotype___AllExpectedGT_pool53/Study_Merge_AllExpectedGT_QYSDCXD23_out.vcf.gz,
    //  /lustre/scratch125/humgen/teams/hgi/mo11/oneK1k/extra0/work/ea/b902eb8a42b41fbeff7b25bd217bc3/Genotype___AllExpectedGT_pool53/Study_Merge_AllExpectedGT_QYSDCXD23_out.vcf.gz.csi,
    //  /lustre/scratch125/humgen/teams/hgi/mo11/oneK1k/extra0/work/0f/954bb3fc06578a6b08173d69999634/cellsnp/cellsnp_pool53/cellSNP.cells.vcf.gz]
    // WARN: Input tuple does not match input set cardinality declared by process `MAIN:YASCP:main_deconvolution:match_genotypes:CONCORDANCE_CALCLULATIONS` -- offending value: [pool75, /lustre/scratch125/humgen/teams/hgi/mo11/oneK1k/extra0/work/bc/e68af9af82417f9ca04323774ca9ff/Genotype___GTMatchedSubset_pool75/Study_Merge_GTMatchedSubset_C7UBEWFX8_out.vcf.gz, /lustre/scratch125/humgen/teams/hgi/mo11/oneK1k/extra0/work/bc/e68af9af82417f9ca04323774ca9ff/Genotype___GTMatchedSubset_pool75/Study_Merge_GTMatchedSubset_C7UBEWFX8_out.vcf.gz.csi, /lustre/scratch125/humgen/teams/hgi/mo11/oneK1k/extra0/work/a3/38e972a5f46b017e4e297402ea076c/Genotype___AllExpectedGT_pool75/Study_Merge_AllExpectedGT_KLGYS5TCV_out.vcf.gz, /lustre/scratch125/humgen/teams/hgi/mo11/oneK1k/extra0/work/a3/38e972a5f46b017e4e297402ea076c/Genotype___AllExpectedGT_pool75/Study_Merge_AllExpectedGT_KLGYS5TCV_out.vcf.gz.csi, /lustre/scratch125/humgen/teams/hgi/mo11/oneK1k/extra0/work/19/1e26bae806d7baf3bf4998fac5f8db/stats_pool75_gt_donor_assignments.csv

    output:
        path("cell_concordance_table.tsv", emit: concordances)

    script:

        """
            echo ${pool_id}
            concordance_calculations_donor_exclusive.py --cpus $task.cpus --cell_vcf ${cell_vcf} --donor_assignments ${donor_table} --gt_match_vcf ${vcf_gt_match} --expected_vcf ${vcf_exp} --cell_assignments ${cell_assignments}
        """
}
