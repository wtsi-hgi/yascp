process gather_handover_dataset {

  publishDir  path: "${outdir}",
              mode: "${params.copy_mode}",
              overwrite: "true"

  when:
    params.data_handover.run_process

  input:
    path(execution_trigger)
    val(outdir_prev)
    path(file__anndata_merged)
    path(file__cellranger_raw_files_table_tsv)
    path(file__cellbender_files_table_tsv)
    path(file__deconv_files_table_tsv)
    path(multiplet_calls)
    path(deconvolution_path)
    path(qc_output_dir)

  output:
    path("${subdir}/*", emit:outfiles_dataset)
    path("${subdir2}/*", emit:outfiles_dataset2)
    path("${subdir2}/adata_celltypes.h5ad", emit:adata_celltypes)
    val(outdir, emit: outdir_dataset)

  script:
    outdir = "${outdir_prev}/handover"
    subdir = "minimal_dataset"
    subdir2="${subdir}_summary"
    if (multiplet_calls) {
      argstr = " --scrublet-output-dir=${qc_output_dir}/multiplet.method=scrublet"
    } else {
      argstr = ""
    }
    """
      gather_minimal_dataset.py \
        --output-dir=${subdir} \
        --cellranger-rawfiles-table=${file__cellranger_raw_files_table_tsv} \
        --cellbender-files-table=${file__cellbender_files_table_tsv} \
        --deconvolution-files-table=${file__deconv_files_table_tsv} \
        --deconvolution-output-dir=${deconvolution_path}/results/split_donor_h5ad \
        --azimuth-output-dir=${qc_output_dir}/azimuth \
        --qc-merged-h5ad=${file__anndata_merged} ${argstr} \
    """
}
