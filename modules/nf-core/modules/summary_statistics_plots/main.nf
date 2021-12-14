process SUMMARY_STATISTICS_PLOTS {
        
    label 'process_low'
    publishDir "${params.outdir}/deconvolution/split_donor_h5ad/${sample}/", mode: "${params.copy_mode}", overwrite: true,
	  saveAs: {filename -> filename.replaceFirst("outputs/","") }
    

    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "/software/hgi/containers/scrna_deconvolution_latest.img"
    } else {
        container "quay.io/biocontainers/multiqc:1.10.1--py_0"
    }

    input: 
      tuple val(sample), val(donor_ids_tsv), path(filtered_matrix_h5), path(scrublet)

    output: 
        
    
    script:
    dir = workflow.workDir+'/../results/deconvolution'
    """

        mkdir -p outputs

        split_h5ad_per_donor.py \\
        --vireo_donor_ids_tsv ${donor_ids_tsv} \\
        --filtered_matrix_h5 ${filtered_matrix_h5} \\
        --scrublet ${scrublet} \\
        --samplename ${sample} \\
        --output_dir \$PWD/outputs \\
        --input_h5_genome_version ${params.split_h5ad_per_donor.input_h5_genome_version} \\
        --print_modules_version ${params.split_h5ad_per_donor.print_modules_version} \\
        --plot_n_cells_per_vireo_donor ${params.split_h5ad_per_donor.plot_n_cells_per_vireo_donor} \\
        --write_donor_level_filtered_cells_h5 ${params.split_h5ad_per_donor.write_donor_level_filtered_cells_h5} \\
        --plotnine_dpi ${params.split_h5ad_per_donor.plotnine_dpi} \\
        --anndata_compression_level ${params.split_h5ad_per_donor.anndata_compression_level}

        # sample h5ad filepath to tsv:
        printf \"$sample\\t\$(find outputs -maxdepth 1 -name '*.h5ad')\" > ${sample}.h5ad.tsv
  
        sed -i 's|outputs/|$dir/split_donor_h5ad/${sample}/|g' ${sample}.h5ad.tsv 

        # deconvoluted donors h5ad file paths to tsv:
        find outputs/donor_level_anndata -maxdepth 1 -name '*.h5ad' -type f -printf \"%f\\n\" | sort | cut -f1 -d'.' > donors.list
        find outputs/donor_level_anndata -maxdepth 1 -name '*.h5ad' | sort > donors.h5ad.list
        paste donors.list donors.h5ad.list > ${sample}.donors.h5ad.tsv
        # paste sample and donor columns 1 and 2 with __
        sed s\"/^/${sample}__/\"g ${sample}.donors.h5ad.tsv > ${sample}__donors.h5ad.tsv
        sed -i 's|outputs/|$dir/split_donor_h5ad/${sample}/|g' ${sample}__donors.h5ad.tsv
        sed -i s\"/^/$sample\\t/\"g ${sample}.donors.h5ad.tsv 
        sed -i 's|outputs/|$dir/split_donor_h5ad/${sample}/|g' ${sample}.donors.h5ad.tsv
        #rm donors.list
        #rm donors.h5ad.list

        # ignore unassigned/doublet h5ad (i.e. assigned cells only):
        cat ${sample}.donors.h5ad.tsv | grep -v unassigned | grep -v doublet > ${sample}.donors.h5ad.assigned.tsv

        cat ${sample}__donors.h5ad.tsv | grep -v unassigned | grep -v doublet > ${sample}__donors.h5ad.assigned.tsv
    """
}
