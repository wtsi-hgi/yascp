process KERAS_CELLTYPE {

    tag { "${experiment_id}" }
    publishDir  path: "${params.outdir}/celltype_assignemt/keras_celltype/${experiment_id}/",
        mode: "${params.copy_mode}",
        overwrite: "true"

    label 'process_high_memory'
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "${params.yascp_container}"

    } else {
        container "${params.yascp_container_docker}"
    }
    input:
        tuple(
            val(experiment_id),
            path(keras_input_h5ad) // from cellbender into scrublet
        )
        path(keras_model)
        path(keras_weights_df)
    
    output:
        tuple(
            val(experiment_id),
            path("*.png"),
            path("*_gtr_*.h5ad"),// do not capture input Pilot_study_of_dissociation_methods_for_human_gut_tissues7980357-scrublet.h5ad
            path("*predictions.h5ad"), 
            emit: out1)
        tuple(
            val(experiment_id),
            path("*predictions.h5ad"), 
            emit: to_merge)
        path('*_celltypes.tsv', emit:predicted_celltype_labels)
        path "versions.yml", emit: versions
    
    script:
        """
            0057-predict_clusters_keras_model-anndata.py \\
            --h5_anndata \"${keras_input_h5ad}\" \\
            --h5_layer \"${params.celltype_prediction.keras.h5_layer}\" \\
            --keras_model \"${keras_model}\" \\
            --keras_weights_df \"${keras_weights_df}\" \\
            --keras_model_cluster_labels \"${params.celltype_prediction.keras.keras_model_cluster_labels}\" \\
            --filter_top_cell_probabilities \"${params.celltype_prediction.keras.filter_top_cell_probabilities}\" \\
            ${params.celltype_prediction.keras.save_all_probabilities} \\
            --output_file \"${experiment_id}___cellbender_fpr${params.cellbender_resolution_to_use}-scrublet-ti_freeze003_prediction\" 

            cat <<-END_VERSIONS > versions.yml
            "${task.process}":
                python: \$(python --version | sed 's/Python //g')
                keras: \$(python -c "import keras; print(keras.__version__)")
                scanpy: \$(python -c "import scanpy; print(scanpy.__version__)")
                argparse: \$(python -c "import argparse; print(argparse.__version__)")
                os: \$(python -c "import os; print(os.__version__)")
                random: \$(python -c "import random; print(random.__version__)")
                warnings: \$(python -c "import warnings; print(warnings.__version__)")
                numpy: \$(python -c "import numpy; print(numpy.__version__)")
                scipy: \$(python -c "import scipy; print(scipy.__version__)")
                pandas: \$(python -c "import pandas; print(pandas.__version__)")
                sklearn: \$(python -c "import sklearn; print(sklearn.__version__)")
                plotnine: \$(python -c "import plotnine; print(plotnine.__version__)")
                matplotlib: \$(python -c "import matplotlib; print(matplotlib.__version__)")
            END_VERSIONS
        """
}
