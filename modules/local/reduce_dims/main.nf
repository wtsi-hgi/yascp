
process DONT_INTEGRATE{

    scratch false      // use tmp directory
    label 'process_medium'
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "${params.yascp_container}"       
    } else {
        container "${params.yascp_container_docker}"
    }

    input:
        val(outdir_prev)
        path(file__anndata)
        path(file__metadata)
        path(file__pcs)
        val(params__pcs)
        each n_pcs

    output:
        val(outdir, emit: outdir)
        path("reduced_dims.tsv.gz", emit: reduced_dims)
        path("${outfile}-dontIntegrate.h5ad", emit: anndata)
        // NOTE: passing the param details as an unpublished file is messy,
        // but I could not get collect of ${param_details} and file to work.
        path(
            "reduced_dims-${param_details}.tsv.gz",
            emit: reduced_dims_params
        )
        path "versions.yml", emit: versions

    script:
        param_details = "${params__pcs}-dontIntegrate"
        param_details = "${param_details}.n_pcs=${n_pcs}"
        outdir = "${outdir_prev}/reduced_dims-${param_details}"
        outfile = "outfile_adata"
        """
        
        0045-reduce_dims_file.py \
            --h5_anndata ${file__anndata} \
            --n_pcs ${n_pcs} \
            --output_file ${outfile}-dontIntegrate
        cp ${outfile}-dontIntegrate-reduced_dims.tsv.gz \
            reduced_dims-${param_details}.tsv.gz
        mv ${outfile}-dontIntegrate-reduced_dims.tsv.gz \
            reduced_dims.tsv.gz

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            python: \$(python --version | sed 's/Python //g')
            bbknn: \$(python -c "import bbknn; print(bbknn.__version__)")
            scanpy: \$(python -c "import scanpy; print(scanpy.__version__)")
            argparse: \$(python -c "import argparse; print(argparse.__version__)")
            distutils: \$(python -c "import distutils; print(distutils.__version__)")
            os: \$(python -c "import os; print(os.__version__)")
            random: \$(python -c "import random; print(random.__version__)")
            numpy: \$(python -c "import numpy; print(numpy.__version__)")
            pandas: \$(python -c "import pandas; print(pandas.__version__)")
        END_VERSIONS
        """

}
