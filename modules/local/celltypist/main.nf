
process CELLTYPIST {
    tag "${model}_${sample}"
    label 'process_medium_memory'
    publishDir "${params.outdir}/celltype_assignemt/celltypist/${model}/${sample}/", mode: "${params.celltypist.copy_mode}", overwrite: true,
	  saveAs: {filename -> filename.replaceFirst("outputs/","").replaceFirst("figures/","") }
    
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {        
      container "${params.yascp_container}"
    } else {
        container "${params.yascp_container_docker}"
    }


    input: 
      tuple val(sample), path(filtered_matrix_h5), val(celltypist_model)

    output: 
      tuple val(sample), path("outputs/*_predicted_labels.csv"), emit: sample_predicted_labels_csv
      path("outputs/*_predicted_labels.csv"), emit: predicted_labels
      tuple val(sample), path("outputs/*_probability_matrix.csv"), emit: sample_probability_matrix_csv
      tuple val(sample), path("outputs/*_decision_matrix.csv"), emit: sample_decision_matrix_csv
      tuple val(sample), path("outputs/*_*.pdf"), emit: sample_plots_pdf
      tuple val(sample), path("outputs/plot_prob/*_*.pdf"), emit: sample_plots_prob_pdf, optional: true

    script:
      model="${celltypist_model}".replaceAll(/^.*[\\/]/, "").replaceFirst(".pkl","")

      if (params.celltypist.sample_plot_probs){
        sample_plot_probs = "--sample_plot_probs"
      }
      else{
        sample_plot_probs = ""
      }

      """

        umask 2 # make files group_writable 
        mkdir -p outputs
        run_celltypist.py \\
          --samplename ${sample} \\
          ${sample_plot_probs} \\
          --filtered_matrix_h5 ${filtered_matrix_h5} \\
          --celltypist_model ${celltypist_model}  \\
          --output_dir \$PWD/outputs  \\
          --input_h5_genome_version ${params.split_h5ad_per_donor.input_h5_genome_version}
      """
}
