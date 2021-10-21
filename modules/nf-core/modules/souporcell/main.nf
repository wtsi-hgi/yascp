
process SOUPORCELL {
    tag "${samplename}"
    label 'process_medium'
    publishDir "${params.outdir}/souporcell/",
    mode: "${params.souporcell.copy_mode}",
    pattern: "${samplename}",
    overwrite: true


    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "/software/hgi/containers/shub_wheaton5_souporcell_latest.img"
    } else {
        log.info 'wrong container, please change this'
        container "quay.io/biocontainers/multiqc:1.10.1--py_0"
    }

    when: 
      params.souporcell.run

    input: 
      tuple val(samplename), path(bam_file), path(bai_file), path(barcodes_tsv_gz), val(souporcell_n_clusters), path(genotypes),val(donors)
      file(reference_fasta)
      

    output:
      tuple val(samplename), file("${samplename}"), emit: souporcell_output_dir
      tuple val(samplename), file("${samplename}/clusters.tsv"), file("${samplename}/cluster_genotypes.vcf"), file("${samplename}/ambient_rna.txt"), emit: souporcell_output_files


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

        souporcell_pipeline.py \\
          -i ${bam_file} \\
          -b bar_codes.txt \\
          -f ${reference_fasta} \\
          -t ${task.cpus} \\
          -o ${samplename} \\
          -k ${souporcell_n_clusters} \\
          --known_genotypes ${genotypes} \\
          --known_genotypes_sample_names ${donors}
      """
}
