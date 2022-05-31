process REPLACE_GT_DONOR_ID{


    publishDir  path: "${params.outdir}/deconvolution/vireo_gt_fix/${samplename}/",
          pattern: "GT_replace_*",
          mode: "${params.copy_mode}",
          overwrite: "true"
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "/software/hgi/containers/mercury_scrna_deconvolution_62bd56a-2021-12-15-4d1ec9312485.sif"
        //// container "/software/hgi/containers/mercury_scrna_deconvolution_latest.img"
    } else {
        container "mercury/scrna_deconvolution:62bd56a"
    }

  label 'process_medium'

  input:
    tuple val(samplename), path(gt_donors), path(vireo_sample_summary),path(vireo___exp_sample_summary),path(vireo__donor_ids),path(vcf_file),path(donor_gt_csi)
    path(gt_match_results)
  output:
    path("test.out", emit: replacements)
    tuple val(samplename), path("GT_replace_donor_ids.tsv"), emit: sample_donor_ids
    tuple val(samplename), path("GT_replace_GT_donors.vireo.vcf.gz"), path(vcf_file),path(donor_gt_csi), emit: sample_donor_vcf
    path("GT_replace_${samplename}.sample_summary.txt"), emit: sample_summary_tsv
    path("GT_replace_${samplename}__exp.sample_summary.txt"), emit: sample__exp_summary_tsv
    path("GT_replace_${samplename}_assignments.tsv"), emit: assignments
    
  script:
    if(params.genotype_phenotype_mapping_file==''){
      in=""
    }else if (params.use_phenotype_ids_for_gt_match){
      in="--genotype_phenotype_mapping ${params.genotype_phenotype_mapping_file}"
      // in=""
    }else{
      in=""
    }

    """
      echo ${samplename} > test.out
      gunzip -k -d --force GT_donors.vireo.vcf.gz
      replace_donors.py -id ${samplename} ${in} --input_file ${params.input_data_table}
      bgzip GT_replace_GT_donors.vireo.vcf
    """
}

process REPLACE_GT_ASSIGNMENTS_WITH_PHENOTYPE{
  label 'process_low'
  publishDir  path: "${params.outdir}/gtmatch/",
          pattern: "*_assignments.csv",
          mode: "${params.copy_mode}",
          overwrite: "true"

  if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
      container "/software/hgi/containers/mercury_scrna_deconvolution_62bd56a-2021-12-15-4d1ec9312485.sif"
      //// container "/software/hgi/containers/mercury_scrna_deconvolution_latest.img"
  } else {
      container "mercury/scrna_deconvolution:62bd56a"
  }

  input:
    path(gt_match_results)

  output:
    path(gt_match_results, emit: donor_match_table)

  script:
    """
      perform_replacement.py --genotype_phenotype_mapping ${params.genotype_phenotype_mapping_file} --assignemts ${gt_match_results}
      
    """

}

process ENHANCE_VIREO_METADATA_WITH_DONOR{
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "/software/hgi/containers/mercury_scrna_deconvolution_62bd56a-2021-12-15-4d1ec9312485.sif"
        //// container "/software/hgi/containers/mercury_scrna_deconvolution_latest.img"
    } else {
        container "mercury/scrna_deconvolution:62bd56a"
    }
  label 'process_small'


  input:
    path(extra_sample_metadata)
    path(donor_n_cells)
    path(out_gt)

  output:
    path('replaced_vireo_exp__donor_n_cells_out.tsv'), emit: replaced_vireo_exp__donor_n_cells_out

  script:
    """
      enhance_vireo_with_metadata.py --Extra_Metadata_Donors ${extra_sample_metadata} --vireo_data ${donor_n_cells}
    """
}



process MATCH_GT_VIREO {
  tag "${pool_id}"

  publishDir  path: "${params.outdir}/gtmatch/",
          pattern: "*_assignments.csv",
          mode: "${params.copy_mode}",
          overwrite: "true"

  if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
      // println "container: /software/hgi/containers/wtsihgi-nf_genotype_match-1.0.sif\n"
      container "/software/hgi/containers/wtsihgi-nf_yascp_htstools-1.0.sif"
  } else {
      container "mercury/wtsihgi-nf_yascp_htstools-1.0"
  }
  label 'process_long'
  //when: params.vireo.run_gtmatch_aposteriori

  input:
    tuple val(pool_id), path(vireo_gt_vcf), path(ref_gt_vcf), path(ref_gt_csi)

  output:
    path("${donor_assignment_csv}", emit: donor_match_table)
    path("${gt_check_output_txt}", emit: gtcheck_out)

  script:
    donor_assignment_csv = "${pool_id}_assignments.csv"
    gt_check_output_txt = "${pool_id}_gtcheck.txt"
  """
    # fix header of vireo VCF
    #tabix -p vcf ${ref_gt_vcf}
    bcftools view -h ${vireo_gt_vcf} > ${pool_id}_header.txt
    bcftools view -Ov ${vireo_gt_vcf} > viewed.vcf
    sed -i '/^##fileformat=VCFv.*/a ##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">' ${pool_id}_header.txt
    bcftools reheader -h ${pool_id}_header.txt -o ${pool_id}_GT_donors.vireo.headfix.vcf.gz viewed.vcf

    # sort and index vireo VCF file (bcftools sort bails out with an error)
    bcftools view ${pool_id}_GT_donors.vireo.headfix.vcf.gz | \
      awk '\$1 ~ /^#/ {print \$0;next} {print \$0 | "sort -k1,1V -k2,2n"}' > ${pool_id}_GT_donors.vireo.srt.vcf
    
    awk '{if(\$0 !~ /^#/) print "chr"\$0; else print \$0}' ${pool_id}_GT_donors.vireo.srt.vcf > ${pool_id}_chr_GT_donors.vireo.srt.vcf
    bgzip ${pool_id}_chr_GT_donors.vireo.srt.vcf
    tabix -p vcf ${pool_id}_chr_GT_donors.vireo.srt.vcf.gz

    bcftools gtcheck -g ${ref_gt_vcf} ${pool_id}_chr_GT_donors.vireo.srt.vcf.gz > ${gt_check_output_txt}

    # generate assignment table
    gtcheck_assign.py ${gt_check_output_txt} ${donor_assignment_csv}

  """
}
