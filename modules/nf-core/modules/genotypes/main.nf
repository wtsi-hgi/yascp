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
