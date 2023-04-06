
process SPLIT_CITESEQ_GEX {
    label 'process_medium'
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "https://yascp.cog.sanger.ac.uk/public/singularity_images/wtsihgi_nf_scrna_qc_6bb6af5-2021-12-23-3270149cf265.sif"
    } else {
        container "wtsihgi/nf_scrna_qc:6bb6af5"
    }

    publishDir  path: "${params.outdir}/citeseq/${sample_name}",
      saveAs: {filename ->
        if (filename.equalsIgnoreCase("antibody-${sample_name}.h5ad")) {
            filename
        }else {
            null
        }
      },
      mode: "${params.copy_mode}",
      overwrite: "true"

    input:
        tuple val(sample_name),path(cellranger_raw) 

    output:
        tuple val(sample_name), path("${sample_name}__gex_data"), emit:gex_data
        tuple val(sample_name), path("antibody-${sample_name}.h5ad"), emit: ab_data2 optional true
        tuple val(sample_name), path("${sample_name}__ab_data"), emit: ab_data optional true

    script:

        """
            strip_citeseq.py --raw_data ${cellranger_raw} -o ${sample_name}
        """
}


process DSB {
    label 'process_medium'
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "https://yascp.cog.sanger.ac.uk/public/singularity_images/wtsihgi_nf_scrna_qc_6bb6af5-2021-12-23-3270149cf265.sif"
    } else {
        container "wtsihgi/nf_scrna_qc:6bb6af5"
    }

    // publishDir  path: "${params.outdir}/citeseq/${sample_name}",
    //   saveAs: {filename ->
    //     if (filename.equalsIgnoreCase("antibody-${sample_name}.h5ad")) {
    //         filename
    //     }else {
    //         null
    //     }
    //   },
    //   mode: "${params.copy_mode}",
    //   overwrite: "true"

    input:
        tuple val(sample_name),path(cellbender_filtered), path(antibody_data), path(cellranger_raw)

    output:
        // tuple val(sample_name), path("${sample_name}__gex_data"), emit:gex_data
        path("output.tmp"), emit: ab_data optional true

    script:
        // [LRBA_B1_BM, /lustre/scratch123/hgi/mdt2/teams/hgi/mo11/tmp_projects/ania/analysis_trego/work/c5/b3f8dd3df7e3cc4c22a8dbd8102e17/cellbender_FPR_0.1_filtered.h5, /lustre/scratch123/hgi/mdt2/teams/hgi/mo11/tmp_projects/ania/analysis_trego/work/27/723d08a62a12435a74f99bab3a683e/antibody-LRBA_B1_BM.h5ad, /lustre/scratch123/hgi/mdt2/teams/hgi/mo11/tmp_projects/ania/analysis_trego/work/27/723d08a62a12435a74f99bab3a683e/LRBA_B1_BM__gex_data]
        """
            echo ${sample_name}
            echo ${cellbender_filtered}
            echo ${antibody_data}
            echo ${cellranger_raw}
            echo 'lets process data with DSB background removal' > output.tmp
        """
}
