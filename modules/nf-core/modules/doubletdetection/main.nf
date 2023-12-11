process doubletdetection {

    tag "${experiment_id}"
    label 'process_low'
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "/lustre/scratch123/hgi/teams/hgi/mo11/tmp_projects/jaguar_yascp/nieks_pipeline/sle-niek/2022-03-07-sc-eQTLgen-pipeline/2022-03-17-scEqtlgen-scripts-and-data/sif_files/WG1/demultiplex/WG1-pipeline-QC_wgpipeline.sif"
    } else {
        container "wtsihgi/nf_scrna_qc:6bb6af5"
    }
    
    publishDir  path: "${params.outdir}/multiplet.method=doubletdetection",
                mode: "${params.copy_mode}",
                overwrite: "true"

    input:
        // path(outdir_prev)
        tuple(
            val(experiment_id),
            path(file_10x_barcodes),
            path(file_10x_features),
            path(file_10x_matrix)
        )

    output:
        path("plots/*.pdf") optional true
        path("plots/*.png") optional true

    script:
        
        outdir = "${params.outdir}/multiplet"
        outdir = "${outdir}.method=doubletdetection"
        outfile = "${experiment_id}"

        """
            TMP_DIR=\$(mktemp -d -p \$(pwd))
            ln --physical ${file_10x_barcodes} \$TMP_DIR
            ln --physical ${file_10x_features} \$TMP_DIR
            ln --physical ${file_10x_matrix} \$TMP_DIR
            DoubletDetection.py --tenxdata_dir \$TMP_DIR --n_iterations 100
        """
}