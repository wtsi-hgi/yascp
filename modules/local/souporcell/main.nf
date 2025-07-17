
process SOUPORCELL {
    tag "${samplename}"
    
    label 'process_high'
    
    publishDir "${params.outdir}/deconvolution/souporcell/",
    mode: "${params.souporcell.copy_mode}",
    pattern: "${samplename}",
    overwrite: true


    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "${params.yascp_container}"
    } else {
        container "${params.yascp_container_docker}"
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

      if (params.genotype_input.run_with_genotype_input){
        known_genotypes = "--known_genotypes ${genotypes}"
        known_genotypes_sample_names ="--known_genotypes_sample_names ${donors}"
      }else{
        known_genotypes = ''
        known_genotypes_sample_names =''
      }

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
          --skip_remap
          -f ${reference_fasta} \\
          -t ${task.cpus} \\
          -o ${samplename} \\
          -k ${souporcell_n_clusters} \\
          ${known_genotypes} \\
          ${known_genotypes_sample_names}
      """
}
