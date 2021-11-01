
process CELLSNP {
    tag "${samplename}"
    
    label 'process_high'
    
    publishDir "${params.outdir}/cellsnp/", mode: "${params.cellsnp.copy_mode}", pattern: "cellsnp_${samplename}", overwrite: true

    
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "/software/hgi/containers/scrna_deconvolution_latest.img"
    } else {
        container "quay.io/biocontainers/multiqc:1.10.1--py_0"
    }

    when: 
        params.cellsnp.run

    input: 
        tuple val(samplename), path(bam_file), path(bai_file), path(barcodes_tsv_gz)
        file(region_vcf)

    output:
    tuple val(samplename), file("cellsnp_${samplename}"), emit: cellsnp_output_dir

    script:
    """
      umask 2 # make files group_writable

      if [[ ${barcodes_tsv_gz} =~ \\.gz\$ ]]; then
        echo \"${barcodes_tsv_gz} is gzipped\"
        zcat ${barcodes_tsv_gz} > bar_codes.txt
      else
        echo \"${barcodes_tsv_gz} is not gzipped\"
        ln -s ${barcodes_tsv_gz} bar_codes.txt
      fi

      cellsnp-lite -s ${bam_file} \\
        -b bar_codes.txt \\
        -O cellsnp_${samplename} \\
        -R ${region_vcf} \\
        -p ${params.cellsnp.p} \\
        --minMAF ${params.cellsnp.min_maf} \\
        --minCOUNT ${params.cellsnp.min_count} --gzip
    """
}
