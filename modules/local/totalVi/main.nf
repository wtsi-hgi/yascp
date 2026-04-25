process TOTAL_VI_INTEGRATION{
    
    if (params.utilise_gpu){
        label 'process_low'
    }else{
        label 'process_medium'
    }
    memory { 
            sizeInGB = adata.size() / 1e9 * 2 * task.attempt
            return (sizeInGB ).toString() + 'GB' 
        }

    publishDir  path: "${outdir_prev}/totalVi",
                saveAs: { filename -> 
                    filename == 'versions.yml' ? null : filename 
                },
                mode: "${params.copy_mode}",
                overwrite: "true"

    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "${params.yascp_container}"
    } else {
        container "${params.yascp_container_docker}"
    }

    input:
        path(adata)
        path(citedata)
        val(outdir_prev)

    output:
        path("./figures"), emit: figs, optional: true
        path("./scvi_model"), emit: scvi_model, optional: true
        path("./totalVI_integrated.h5ad"), emit: totalVI_integrated, optional: true
        path "versions.yml", emit: versions

        
    script:

        """
            totalVI.py -h5ad_file ${adata}
            cat <<-END_VERSIONS > versions.yml
            "${task.process}":
                python: \$(python --version | sed 's/Python //g')
                python library argparse: \$(python -c "import argparse; print(argparse.__version__)")
                python library matplotlib: \$(python -c "import matplotlib; print(matplotlib.__version__)")
                python library numpy: \$(python -c "import numpy; print(numpy.__version__)")
                python library pandas: \$(python -c "import pandas; print(pandas.__version__)")
                python library scanpy: \$(python -c "import scanpy; print(scanpy.__version__)")
                python library scipy: \$(python -c "import scipy; print(scipy.__version__)")
                python library scvi: \$(python -c "import scvi; print(scvi.__version__)")
            END_VERSIONS
        """

}
