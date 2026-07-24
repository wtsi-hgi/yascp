process CONVERT_H5AD_TO_MTX {
    tag "${samplename}"    
    label 'process_low'
   
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "${params.yascp_container}"
    } else {
        container "wtsihgi/nf_scrna_qc:6bb6af5"
    }

    input:
        path(h5ad)
    output:
        tuple val("${h5ad.getBaseName()}"), path("${h5ad.getBaseName()}"), emit: channel__file_paths_10x
        path "versions.yml", emit: versions

    script:
        """
        # Remove the .h5ad extension to get the base name
        base_name=\$(basename ${h5ad} .h5ad)
        
        # Create the output directory with the base name
        out_dir=\${base_name}
        mkdir -p \$out_dir

        # Run the Python script with the base name as the out_file prefix
        h5ad_to_tenxmatrix.py ${h5ad} \$out_dir
        
        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            python: \$(python --version | sed 's/Python //g')
            python library anndata: \$(python -c "import anndata; print(anndata.__version__)")
            python library argparse: \$(python -c "import argparse; print(argparse.__version__)")
            python library distutils: \$(python -c "import distutils; print(distutils.__version__)")
            python library numpy: \$(python -c "import numpy; print(numpy.__version__)")
            python library pandas: \$(python -c "import pandas; print(pandas.__version__)")
            python library scanpy: \$(python -c "import scanpy; print(scanpy.__version__)")
            python library scipy: \$(python -c "import scipy; print(scipy.__version__)")
        END_VERSIONS
        
        """
}


process CONVERT_MTX_TO_H5AD {
    tag "${samplename}"    
    label 'process_low'
   
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "${params.yascp_container}"
    } else {
        container "wtsihgi/nf_scrna_qc:6bb6af5"
    }

    input:
        tuple val(name), path(mtx1)
    output:
        tuple val(name), path("${name}.h5ad"), emit: gex_h5ad
        path "versions.yml", emit: versions

    script:
        """
        # Run the Python script to convert MTX to H5AD
        h5ad_from_tenxmatrix.py ${mtx1} ${name}.h5ad
        
        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            python: \$(python --version | sed 's/Python //g')
            python library anndata: \$(python -c "import anndata; print(anndata.__version__)")
            python library argparse: \$(python -c "import argparse; print(argparse.__version__)")
            python library pandas: \$(python -c "import pandas; print(pandas.__version__)")
            python library scanpy: \$(python -c "import scanpy; print(scanpy.__version__)")
            python library scipy: \$(python -c "import scipy; print(scipy.__version__)")
        END_VERSIONS


        """
}