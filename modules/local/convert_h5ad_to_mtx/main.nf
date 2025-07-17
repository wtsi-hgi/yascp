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

    script:
        """
        # Remove the .h5ad extension to get the base name
        base_name=\$(basename ${h5ad} .h5ad)
        
        # Create the output directory with the base name
        out_dir=\${base_name}
        mkdir -p \$out_dir

        # Run the Python script with the base name as the out_file prefix
        h5ad_to_tenxmatrix.py ${h5ad} \$out_dir
        
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

    script:
        """
        # Run the Python script to convert MTX to H5AD
        h5ad_from_tenxmatrix.py ${mtx1} ${name}.h5ad
        """
}