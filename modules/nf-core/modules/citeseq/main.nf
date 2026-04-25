
process SPLIT_CITESEQ_GEX {
    label 'process_medium'
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "https://yascp.cog.sanger.ac.uk/public/singularity_images/wtsihgi_nf_scrna_qc_6bb6af5-2021-12-23-3270149cf265.sif"
    } else {
        container "wtsihgi/nf_scrna_qc:6bb6af5"
    }

    publishDir  path: "${params.outdir}/citeseq/${sample_name}",
      saveAs: {filename ->
        if (filename.contains("antibody-")) {
            filename.replaceAll("antibody-", "${mode}_antibody-")
        }else {
            null
        }
      },
      mode: "${params.copy_mode}",
      overwrite: "true"

    input:
        tuple val(sample_name),path(cellranger_raw) 
        val(mode)

    output:
        tuple val(sample_name), path("${sample_name}__gex_data"), emit:gex_data
        tuple val(sample_name), path("antibody-${sample_name}.h5ad"), emit: ab_data2 optional true
        tuple val(sample_name), path("${sample_name}__ab_data"), emit: ab_data optional true
        tuple val(sample_name), path("${sample_name}__gex_data/barcodes.tsv.gz"), path("${sample_name}__gex_data/features.tsv.gz"), path("${sample_name}__gex_data/matrix.mtx.gz"), emit: channel__file_paths_10x
 
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
