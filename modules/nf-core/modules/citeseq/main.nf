
process SPLIT_CITESEQ_GEX {
    label 'process_medium'
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "https://yascp.cog.sanger.ac.uk/public/singularity_images/wtsihgi_nf_scrna_qc_6bb6af5-2021-12-23-3270149cf265.sif"
    } else {
        container "wtsihgi/nf_scrna_qc:6bb6af5"
    }

    publishDir  path: "${params.outdir}/citeseq/${mode}/${sample_name}",
    //   saveAs: {filename ->
    //     if (filename.contains("antibody-")) {
    //         filename.replaceAll("antibody-", "${mode}_antibody-")
    //     }else {
    //         null
    //     }
    //   },
      mode: "${params.copy_mode}",
      overwrite: "true"

    input:
        tuple val(sample_name),path(cellranger_raw) 
        val(mode)

    output:
        tuple val(sample_name), path("${sample_name}__Gene_Expression"), emit:gex_data
        tuple val(sample_name), path("antibody-${sample_name}.h5ad"), emit: ab_data2 optional true
        tuple val(sample_name), path("${sample_name}__*"), emit: ab_data
        tuple val(sample_name), path("${sample_name}__Gene_Expression/barcodes.tsv.gz"), path("${sample_name}__Gene_Expression/features.tsv.gz"), path("${sample_name}__Gene_Expression/matrix.mtx.gz"), emit: channel__file_paths_10x
 
    script:

        """
            strip_citeseq.py --raw_data ${cellranger_raw} -o ${sample_name}
        """
}


process DSB {
    label 'process_medium'
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "https://yascp.cog.sanger.ac.uk/public/singularity_images/azimuth_dsb.img"
    } else {
        container "mercury/azimuth_dsb:latest"
    }

    publishDir  path: "${params.outdir}/citeseq/${sample_name}", mode: "${params.copy_mode}",
      overwrite: "true"

    input:
        tuple val(sample_name),path(cellbender_filtered), path(antibody_data), path(cellranger_raw)


//         filtered_cellranger = sample_filtered_feature_bc_matrix.h5
//         cellranger_rawfile_path = 'raw_feature_bc_matrix.h5'
//         sample_name <-'cellranger700_multi_bc45a1c2fe2a3fbbcde46cf984cf42e2'
//         file_with_qc_applied <-'cellranger700_multi_bc45a1c2fe2a3fbbcde46cf984cf42e2___sample_QCd_adata.h5ad'


    output:
        // tuple val(sample_name), path("${sample_name}__gex_data"), emit:gex_data
        path("*.pdf"), emit: plots optional true
        path("*.dsb_technical_stats.RDS"), emit: dsb_technical_stats optional true
        path("*.dsb_protein_stats.RDS"), emit: dsb_protein_stats optional true
        path("*_firstrun_dsb.h5Seurat"), emit: firstrun_dsb optional true
        path(antibody_data), emit: antibody_data optional true
        path(cellranger_raw), emit: cellranger_raw optional true

    script:
        // [LRBA_B1_BM, /lustre/scratch123/hgi/mdt2/teams/hgi/mo11/tmp_projects/ania/analysis_trego/work/c5/b3f8dd3df7e3cc4c22a8dbd8102e17/cellbender_FPR_0.1_filtered.h5, /lustre/scratch123/hgi/mdt2/teams/hgi/mo11/tmp_projects/ania/analysis_trego/work/27/723d08a62a12435a74f99bab3a683e/antibody-LRBA_B1_BM.h5ad, /lustre/scratch123/hgi/mdt2/teams/hgi/mo11/tmp_projects/ania/analysis_trego/work/27/723d08a62a12435a74f99bab3a683e/LRBA_B1_BM__gex_data]
        """
            echo ${sample_name}
            echo ${cellbender_filtered}
            echo ${antibody_data}
            echo ${cellranger_raw}

            dsb_process_extr.R ${cellbender_filtered} ${cellranger_raw} ${antibody_data} ${sample_name}
            echo 'lets process data with DSB background removal' > output.tmp
        """
}

process DSB_INTEGRATE{

    label 'process_medium'
    tag "${sample_name}"
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "/lustre/scratch123/hgi/teams/hgi/mo11/tmp_projects/jaguar_yascp/nieks_pipeline/yascp_run/azimuth_dsb_17_02_2024.sif"
    } else {
        container "mercury/azimuth_dsb:latest"
    }

    publishDir  path: "${params.outdir}/citeseq/all_data_integrated", mode: "${params.copy_mode}",
      overwrite: "true"
    
    output:
        path("out"), emit: all_data_integrated
        path("./out/tmp_rds_files/all_samples_integrated.RDS"), emit: all_data_integrated_rds

    input:
        path(citeseq_rsd)
        path(vireo)
        path(assignments_all_pools)
        path(tmp_rsd)
        path(matched_donors)

    script:
    """
    echo 'running1'
    2.integrate.R
    """

}

process MULTIMODAL_INTEGRATION{

    label 'process_medium'
    tag "${sample_name}"
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "/lustre/scratch123/hgi/teams/hgi/mo11/tmp_projects/jaguar_yascp/nieks_pipeline/yascp_run/azimuth_dsb_26_02_2024.sif"
    } else {
        container "mercury/azimuth_dsb:latest"
    }

    // publishDir  path: "${params.outdir}/citeseq/all_data_integrated", mode: "${params.copy_mode}",
    //   overwrite: "true"
    
    // output:
    //     path("out"), emit: all_data_integrated
    //     path("./out/tmp_rds_files/all_samples_integrated.RDS"), emit: all_data_integrated_rds

    input:
        path(all_samples_integrated)
        // path(vireo)
        // path(assignments_all_pools)
        // path(tmp_rsd)
        // path(matched_donors)

    script:
    """
    echo 'running1'
    3.WNN_integrate_SCT_CITE.R
    """

}


process DSB_PROCESS {
    label 'process_medium'
    tag "${sample_name}"
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "/lustre/scratch123/hgi/teams/hgi/mo11/tmp_projects/jaguar_yascp/nieks_pipeline/yascp_run/work/ef/1cda7f07670de6124e43c7c63adfea/azimuth_dsb_14_02_2024.img"
    } else {
        container "mercury/azimuth_dsb:latest"
    }

    publishDir  path: "${params.outdir}/citeseq/${sample_name}",      
        saveAs: {filename ->
        if (filename.contains("tmp_rds_files__")) {
            null
        }else {
            filename
        }
      }, mode: "${params.copy_mode}",
      overwrite: "true"

    input:
        tuple val(sample_name),path(cellranger_rawfile_path), path(filtered_feature_bc_matrix), path(sample_QCd_adata)
    


    output:
        path("CITE__*"), emit: citeseq_rsd
        path('tmp_rds_files__*'), emit: tmp_rsd
    script:
        """
  
            add_adt.R ${sample_name} ${cellranger_rawfile_path} ${filtered_feature_bc_matrix} ${sample_QCd_adata}
        """
}