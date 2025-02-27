process CONVERT_H5AD_TO_MTX {
    tag "${samplename}"    
    label 'process_low'
   
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "https://yascp.cog.sanger.ac.uk/public/singularity_images/wtsihgi_nf_scrna_qc_6bb6af5-2021-12-23-3270149cf265.sif"
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
        container "https://yascp.cog.sanger.ac.uk/public/singularity_images/wtsihgi_nf_scrna_qc_6bb6af5-2021-12-23-3270149cf265.sif"
    } else {
        container "wtsihgi/nf_scrna_qc:6bb6af5"
    }

    input:
        tuple val(name), path(mtx1),path(mtx2),path(mtx3)
    output:
        tuple val(name), path("${name}.h5ad"), emit: gex_h5ad

    script:
        """

       mkdir -p work_dir

        # Link the MTX folders to the working directory
        # Extract the base names of the MTX files

        # Link the MTX files to the working directory with the same name
        cd work_dir && ln -s ../${mtx1} && ln -s ../${mtx2} && ln -s ../${mtx3} && cd ..

        # Run the Python script to convert MTX to H5AD
        h5ad_from_tenxmatrix.py work_dir ${name}.h5ad
        """
}