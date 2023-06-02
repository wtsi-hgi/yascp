process capture_cellsnp_files{
  publishDir  path: "${params.outdir}"
  // cache false
  label 'process_tiny'
  input:
    path(cellsnp_location)
   
  output:
    path("output_cellsnp.csv"),emit:cellsnp_loc optional true
    path(cellsnp_location) optional true
  script:
  """
    echo '${params.cellsnp_recapture}'
    for OUTPUT in \$(ls ${cellsnp_location})
    do
    samplename1=\$(echo \$OUTPUT | sed 's/cellsnp_//g') 
    echo "\$samplename1 \$PWD/${cellsnp_location}/\$OUTPUT" >> output_cellsnp.csv
    done
  """    

}


process DYNAMIC_DONOR_EXCLUSIVE_SNP_SELECTION{
    label 'process_medium'
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "https://yascp.cog.sanger.ac.uk/public/singularity_images/mercury_scrna_deconvolution_62bd56a-2021-12-15-4d1ec9312485.sif"
        //// container "https://yascp.cog.sanger.ac.uk/public/singularity_images/mercury_scrna_deconvolution_latest.img"
    } else {
        container "mercury/scrna_deconvolution:62bd56a"
    }

    input: 
        tuple val(samplename), path(vcf_file),path(csi),path(cellsnp_primary_file)
    output:
      tuple val(samplename), path("cellsnp_panel_${samplename}.vcf.gz"),emit:cellsnp_pool_panel
      tuple val(samplename), path("set2_informative_sites_${samplename}.tsv"), path("set2_informative_sites_${samplename}.tsv"),emit:informative_uninformative_sites 
    script:       
      """
        echo ${samplename}
        echo ${vcf_file}
        echo ${cellsnp_primary_file}
        bcftools view -i 'MAF > 0.0001 & R2>=1.00' -Oz -o dynamic_snps.vcf.gz ${vcf_file}
        dynamic_donor_exclusive_snp_selection.py -cpus ${task.cpus} -vcf dynamic_snps.vcf.gz -cellsnp ${cellsnp_primary_file}
        echo test > output.csv
        bcftools view -h ${cellsnp_primary_file} > cellsnp_panel_${samplename}.vcf
        cat cellsnp_variants.tsv >> cellsnp_panel_${samplename}.vcf
        ln -s set1_uninformative_sites.tsv set1_uninformative_sites_${samplename}.tsv
        ln -s set2_informative_sites.tsv set2_informative_sites_${samplename}.tsv
        gzip cellsnp_panel_${samplename}.vcf
        rm -r dynamic_snps.vcf.gz
      """
}



process CELLSNP {
    tag "${samplename}"
    
    label 'many_cores_small_mem'
    
    publishDir "${params.outdir}/cellsnp/", mode: "${params.copy_mode}", pattern: "cellsnp_${samplename}", overwrite: true

    
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "https://yascp.cog.sanger.ac.uk/public/singularity_images/mercury_scrna_deconvolution_62bd56a-2021-12-15-4d1ec9312485.sif"
        //// container "https://yascp.cog.sanger.ac.uk/public/singularity_images/mercury_scrna_deconvolution_latest.img"
    } else {
        container "mercury/scrna_deconvolution:62bd56a"
    }

    input: 
        tuple val(samplename), path(bam_file), path(bai_file), path(barcodes_tsv_gz),val(n_pooled),path(region_vcf)

    output:
      tuple val(samplename), file("cellsnp_${samplename}"), emit: cellsnp_output_dir
      tuple val(samplename), path("cellsnp_${samplename}/cellSNP.cells.vcf.gz"), emit: cell_vcfs

    script:
    if (n_pooled=='1'){
      genotype_file=' --genotype '
    }else{
      genotype_file=' --genotype '
    }
    """
      echo ${n_pooled}
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
        -p ${task.cpus} \\
        --minMAF ${params.cellsnp.min_maf} \\
        --minCOUNT ${params.cellsnp.min_count} --gzip ${genotype_file}
    """
}
