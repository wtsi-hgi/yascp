
process SPLIT_CITESEQ_GEX {
    label 'process_medium'
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "${params.yascp_container}"
    } else {
        container "${params.yascp_container_docker}"
    }

    publishDir  path: "${params.outdir}/preprocessing/data_modalities_split/${mode}/${sample_name}",
        saveAs: { filename -> 
        (filename == 'versions.yml' || filename.endsWith('_barcodes.counts.txt')) ? null : filename 
        },
    mode: "${params.copy_mode}",
    overwrite: "true"

    input:
        tuple val(sample_name),path(cellranger_raw) 
        val(mode)

    output:
        tuple val(sample_name), path("${sample_name}__*"), path("*__Multiplexing_Capture.tsv"), emit: multiplexing_capture_channel_for_demultiplexing  optional true
        tuple val(sample_name), path("${sample_name}__Gene_Expression"), emit:gex_data
        tuple val(sample_name), path("antibody-${sample_name}.h5ad"), emit: ab_data2 optional true
        tuple val(sample_name), path("Gene_Expression-${sample_name}.h5ad"), emit: gex_h5ad optional true
        path("Gene_Expression-${sample_name}.h5ad"), emit: gex_h5ad_2 optional true
        path("*.tsv"), emit: quants_data optional true
        tuple val(sample_name), path("${sample_name}__*"), emit: ab_data optional true
        tuple val(sample_name), path("${sample_name}__Gene_Expression/barcodes.tsv.gz"), path("${sample_name}__Gene_Expression/features.tsv.gz"), path("${sample_name}__Gene_Expression/matrix.mtx.gz"), emit: channel__file_paths_10x
        tuple val(sample_name), path("${sample_name}_barcodes.counts.txt"), emit:output_check
        path "versions.yml", emit: versions
 
    script:

        """
            matrix_file="${cellranger_raw}/matrix.mtx"
            compressed_file="${cellranger_raw}/matrix.mtx.gz"

            # Check if the compressed file exists
            if [ ! -f "\$compressed_file" ]; then
                echo "\$compressed_file does not exist. Compressing \$matrix_file..."
            
                # Compress the file without deleting the original
                gzip -c "\$matrix_file" > "\$compressed_file"
            
                echo "Compression complete. \$matrix_file has been compressed to \$compressed_file"
            else
                echo "\$compressed_file already exists. No action needed."
            fi

            matrix_file="${cellranger_raw}/barcodes.tsv"
            compressed_file="${cellranger_raw}/barcodes.tsv.gz"
            # Check if the compressed file exists
            if [ ! -f "\$compressed_file" ]; then
                echo "\$compressed_file does not exist. Compressing \$matrix_file..."
            
                # Compress the file without deleting the original
                gzip -c "\$matrix_file" > "\$compressed_file"
            
                echo "Compression complete. \$matrix_file has been compressed to \$compressed_file"
            else
                echo "\$compressed_file already exists. No action needed."
            fi

            features_file="${cellranger_raw}/features.tsv.gz"
            peaks_file="${cellranger_raw}/peaks.bed"

            # Check if the features.tsv.gz file exists
            if [ ! -f "\$features_file" ]; then
            echo "\$features_file does not exist. Creating it from \$peaks_file..."

            # Create the features.tsv file from the peaks.bed file
            awk 'BEGIN{OFS="\t"} {print \$1 ":" \$2 "-" \$3, \$1 ":" \$2 "-" \$3, "Gene Expression"}' "\$peaks_file" | gzip > "\$features_file"

            echo "Creation of \$features_file is complete."
            else
            echo "\$features_file already exists. No action needed."
            fi

            strip_citeseq.py --raw_data ${cellranger_raw} -o ${sample_name} -ha ${params.citeseq_config.hastag_multiplexing_capture_labels}
            zcat ${sample_name}__Gene_Expression/barcodes.tsv.gz | wc -l > ${sample_name}_barcodes.counts.txt

            cat <<-END_VERSIONS > versions.yml
            "${task.process}":
                python: \$(python --version | sed 's/Python //g')
                python library anndata: \$(python -c "import anndata; print(anndata.__version__)")
                python library argparse: \$(python -c "import argparse; print(argparse.__version__)")
                python library click: \$(python -c "import click; print(click.__version__)")
                python library distutils: \$(python -c "import distutils; print(distutils.__version__)")
                python library logging: \$(python -c "import logging; print(logging.__version__)")
                python library numpy: \$(python -c "import numpy; print(numpy.__version__)")
                python library pandas: \$(python -c "import pandas; print(pandas.__version__)")
                python library scanpy: \$(python -c "import scanpy; print(scanpy.__version__)")
                python library scipy: \$(python -c "import scipy; print(scipy.__version__)")
                python library tables: \$(python -c "import tables; print(tables.__version__)")
            END_VERSIONS

        """
}


process HASTAG_DEMULTIPLEX {
    label 'process_medium'
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "${params.yascp_container}"
    } else {
        container "${params.yascp_container_docker}"
    }

    publishDir  path: "${params.outdir}/deconvolution/hastag_demultiplex/${sample_name}",
    saveAs: { filename -> 
            filename == 'versions.yml' ? null : filename 
        }, 
        mode: "${params.copy_mode}",
      overwrite: "true"

    input:
        tuple val(sample_name), path(paths), path(multiplexing_capture)
    output:
        path("${sample_name}__hastag_demux_results.tsv"), emit: results
        path "versions.yml", emit: versions
    script:
        """
            hastag_demultiplex.R
            ln -s hastag_demux_results.tsv ${sample_name}__hastag_demux_results.tsv
            
            cat <<-END_VERSIONS > versions.yml
            "${task.process}":
                r-base: \$(R --version | sed -n '1p' | sed 's/R version //; s/ (.*//')
                r library Seurat: \$(Rscript -e "cat(as.character(packageVersion('Seurat')))")
                r library tidyverse: \$(Rscript -e "cat(as.character(packageVersion('tidyverse')))")
            END_VERSIONS

        """
}

process DSB_INTEGRATE{

    label 'process_medium'
    tag "${sample_name}"
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "${params.yascp_container}"
    } else {
        container "${params.yascp_container_docker}"
    }

    publishDir  path: "${params.outdir}/citeseq/all_data_integrated",
    saveAs: { filename -> 
            filename == 'versions.yml' ? null : filename 
        }, 
        mode: "${params.copy_mode}",
      overwrite: "true"
    
    output:
        // path("out"), emit: outdir
        path("*all_samples_integrated.RDS"), emit: tmp_rds_file
        path("figures__*/*"), emit: figs
        path "versions.yml", emit: versions
        


    input:
        path(tmp_rsd)
        each vars_to_regress
        val(k_anchor)
        val(dims)
        val(ndim_sct)
        val(ndim_citeBgRemoved)
        val(ndim_cite_integrated)
    script:

        if (vars_to_regress == ''){
            vars_to_regress='NONE'
        }
        
        """
        echo 'running1'
        integrate.R ${vars_to_regress} ${k_anchor} ${dims} ${ndim_sct} ${ndim_citeBgRemoved} ${ndim_cite_integrated}

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            r-base: \$(R --version | sed -n '1p' | sed 's/R version //; s/ (.*//')
            r library future.apply: \$(Rscript -e "cat(as.character(packageVersion('future.apply')))")
            r library ggplot2: \$(Rscript -e "cat(as.character(packageVersion('ggplot2')))")
            r library Matrix: \$(Rscript -e "cat(as.character(packageVersion('Matrix')))")
            r library RColorBrewer: \$(Rscript -e "cat(as.character(packageVersion('RColorBrewer')))")
            r library Seurat: \$(Rscript -e "cat(as.character(packageVersion('Seurat')))")
        END_VERSIONS

        """

}

process MULTIMODAL_INTEGRATION{

    label 'process_medium'
    tag "${sample_name}"
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "${params.yascp_container}"
    } else {
        container "${params.yascp_container_docker}"
    }

    publishDir  path: "${params.outdir}/citeseq/all_data_integrated",
    saveAs: { filename -> 
            filename == 'versions.yml' ? null : filename 
        }, 
        mode: "${params.copy_mode}",
      overwrite: "true"
    
    output:
        path("*wnn.integrated.RDS"), emit: wnn_integrated_file
        path "versions.yml", emit: versions

    input:
        path(tmp_rds_file)

    script:
    """
    echo 'running1'
    WNN_integrate_SCT_CITE.R ${tmp_rds_file}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        r-base: \$(R --version | sed -n '1p' | sed 's/R version //; s/ (.*//')
        r library future.apply: \$(Rscript -e "cat(as.character(packageVersion('future.apply')))")
        r library ggbeeswarm: \$(Rscript -e "cat(as.character(packageVersion('ggbeeswarm')))")
        r library ggplot2: \$(Rscript -e "cat(as.character(packageVersion('ggplot2')))")
        r library patchwork: \$(Rscript -e "cat(as.character(packageVersion('patchwork')))")
        r library RColorBrewer: \$(Rscript -e "cat(as.character(packageVersion('RColorBrewer')))")
        r library Seurat: \$(Rscript -e "cat(as.character(packageVersion('Seurat')))")
        r library tidyverse: \$(Rscript -e "cat(as.character(packageVersion('tidyverse')))")
    END_VERSIONS

    """

}

process VDJ_INTEGRATION{

    label 'process_medium'
    tag "${sample_name}"
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "${params.yascp_container}"
    } else {
        container "${params.yascp_container_docker}"
    }

    publishDir  path: "${params.outdir}/citeseq/all_data_integrated",
    saveAs: { filename -> 
            filename == 'versions.yml' ? null : filename 
        }, 
        mode: "${params.copy_mode}",
      overwrite: "true"
    
    output:
        path("*all_samples_integrated.vdj.RDS"), emit: all_data_integrated_vdj_rds
        path("*all_samples_integrated.BCR.RDS"), emit: all_data_integrated_BCR_rds
        path("*all_samples_integrated.TCR.RDS"), emit: all_data_integrated_TCR_rds
        path "versions.yml", emit: versions

    input:
        path(all_cellranger_samples)
        path(wnn_integrated_file)


    script:
    """
    echo 'running1'
    add_vdj.R ${wnn_integrated_file}
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        r-base: \$(R --version | sed -n '1p' | sed 's/R version //; s/ (.*//')
        r library future.apply: \$(Rscript -e "cat(as.character(packageVersion('future.apply')))")
        r library ggplot2: \$(Rscript -e "cat(as.character(packageVersion('ggplot2')))")
        r library ggpubr: \$(Rscript -e "cat(as.character(packageVersion('ggpubr')))")
        r library immunarch: \$(Rscript -e "cat(as.character(packageVersion('immunarch')))")
        r library RColorBrewer: \$(Rscript -e "cat(as.character(packageVersion('RColorBrewer')))")
        r library Seurat: \$(Rscript -e "cat(as.character(packageVersion('Seurat')))")
        r library tidyverse: \$(Rscript -e "cat(as.character(packageVersion('tidyverse')))")
    END_VERSIONS

    """

}


process PREPROCESS_PROCESS {
    label 'process_medium'
    tag "${sample_name}"
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "${params.yascp_container}"
    } else {
        container "${params.yascp_container_docker}"
    }

    publishDir  path: "${params.outdir}/citeseq/DSB/${sample_name}",      
        saveAs: {filename ->
        if (filename.contains("tmp_rds_files__")) {
            null
        } else if(filename == 'versions.yml') {
            null
        } else {
            filename
        }
      }, mode: "${params.copy_mode}",
      overwrite: "true"

    input:
        tuple val(sample_name), path(vireo_path), path(rds_path),path(matched_donors)
        each vars_to_regress

    output:
        path("normalised__${sample_name}.withADT.RDS"), emit: tmp_rsd
        path "versions.yml", emit: versions

    script:


        if (vars_to_regress == ''){
            vars_to_regress='NONE'
        }
        """
            process_donor_data_for_integration.R ${sample_name} ${vireo_path} ${matched_donors} ${rds_path} ${vars_to_regress}
            cat <<-END_VERSIONS > versions.yml
            "${task.process}":
                r-base: \$(R --version | sed -n '1p' | sed 's/R version //; s/ (.*//')
                r library future.apply: \$(Rscript -e "cat(as.character(packageVersion('future.apply')))")
                r library ggplot2: \$(Rscript -e "cat(as.character(packageVersion('ggplot2')))")
                r library Matrix: \$(Rscript -e "cat(as.character(packageVersion('Matrix')))")
                r library RColorBrewer: \$(Rscript -e "cat(as.character(packageVersion('RColorBrewer')))")
                r library Seurat: \$(Rscript -e "cat(as.character(packageVersion('Seurat')))")
            END_VERSIONS

        """
}


process DSB_PROCESS {
    label 'process_medium'
    tag "${sample_name}"
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "${params.yascp_container}"
    } else {
        container "${params.yascp_container_docker}"
    }

    publishDir  path: "${params.outdir}/citeseq/DSB/${sample_name}",      
        saveAs: {filename ->
        if (filename.contains("tmp_rds_files__")) {
            null
        } else if(filename == 'versions.yml') {
            null
        } else {
            filename
        }
      }, mode: "${params.copy_mode}",
      overwrite: "true"

    input:
        tuple val(sample_name),path(cellranger_rawfile_path), path(filtered_feature_bc_matrix), path(sample_QCd_adata)
    


    output:
        path("CITE__*"), emit: citeseq_rsd
        path("tmp_rds_files__*/*/${sample_name}*.RDS"), emit: tmp_rsd
        tuple val(sample_name), path("tmp_rds_files__*/*/${sample_name}*.RDS"), emit: ch_for_norm
        path "versions.yml", emit: versions
    script:
        """
   
            add_adt.R ${sample_name} ${cellranger_rawfile_path} ${filtered_feature_bc_matrix} ${sample_QCd_adata}
            
            cat <<-END_VERSIONS > versions.yml
            "${task.process}":
                r-base: \$(R --version | sed -n '1p' | sed 's/R version //; s/ (.*//')
                r library dsb: \$(Rscript -e "cat(as.character(packageVersion('dsb')))")
                r library future.apply: \$(Rscript -e "cat(as.character(packageVersion('future.apply')))")
                r library ggplot2: \$(Rscript -e "cat(as.character(packageVersion('ggplot2')))")
                r library ggpubr: \$(Rscript -e "cat(as.character(packageVersion('ggpubr')))")
                r library RColorBrewer: \$(Rscript -e "cat(as.character(packageVersion('RColorBrewer')))")
                r library Seurat: \$(Rscript -e "cat(as.character(packageVersion('Seurat')))")
                r library SeuratDisk: \$(Rscript -e "cat(as.character(packageVersion('SeuratDisk')))")
                r library tidyverse: \$(Rscript -e "cat(as.character(packageVersion('tidyverse')))")
                r library viridis: \$(Rscript -e "cat(as.character(packageVersion('viridis')))")
            END_VERSIONS

        """
}