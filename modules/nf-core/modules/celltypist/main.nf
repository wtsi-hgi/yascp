
process CELLTYPIST {
    tag "${samplename}"
    label 'process_high_memory'
    publishDir "${params.outdir}/celltype/celltypist/${model}/${sample}/", mode: "${params.celltypist.copy_mode}", overwrite: true,
	  saveAs: {filename -> filename.replaceFirst("outputs/","").replaceFirst("figures/","") }
    
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "https://yascp.cog.sanger.ac.uk/public/singularity_images/mercury_scrna_deconvolution_62bd56a-2021-12-15-4d1ec9312485.sif"
        //// container "https://yascp.cog.sanger.ac.uk/public/singularity_images/mercury_scrna_deconvolution_latest.img"
    } else {
        container "mercury/scrna_deconvolution:62bd56a"
    }

    when: 
      params.celltypist.run

    input: 
      tuple val(sample), val(filtered_matrix_h5), val(celltypist_model)

    output: 
      tuple val(sample), path("outputs/*_predicted_labels.csv"), emit: sample_predicted_labels_csv
      path("outputs/*_predicted_labels.csv"), emit: predicted_labels
      tuple val(sample), path("outputs/*_probability_matrix.csv"), emit: sample_probability_matrix_csv
      tuple val(sample), path("outputs/*_decision_matrix.csv"), emit: sample_decision_matrix_csv
      tuple val(sample), path("outputs/*_*.pdf"), emit: sample_plots_pdf
      tuple val(sample), path("outputs/plot_prob/*_*.pdf"), emit: sample_plots_prob_pdf

    script:
      model="${celltypist_model}".replaceFirst(".pkl","")

      filtered_matrix_h5_path = file("${filtered_matrix_h5}/../filtered_feature_bc_matrix.h5")
      if (filtered_matrix_h5_path.exists()){
        _=""
        
      }else{
        filtered_matrix_h5_path = file("${filtered_matrix_h5}/../cellbender_FPR_0pt05_filtered.h5")
      }


      """

        umask 2 # make files group_writable 
        mkdir -p outputs
        run_celltypist.py \\
          --samplename ${sample} \\
          --filtered_matrix_h5 ${filtered_matrix_h5} \\
          --celltypist_model ${celltypist_model}  \\
          --output_dir \$PWD/outputs  \\
          --input_h5_genome_version ${params.split_h5ad_per_donor.input_h5_genome_version}
      """
}
