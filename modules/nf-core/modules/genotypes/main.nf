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

  script:
    """

      echo ${samplename} > test.out
      gunzip -k -d --force GT_donors.vireo.vcf.gz
      replace_donors.py -id ${samplename}
      bgzip GT_replace_GT_donors.vireo.vcf
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
      container "/software/hgi/containers/wtsihgi-nf_yascp_htstools-1.1.sif"
  } else {
      container "mercury/wtsihgi-nf_yascp_htstools-1.1"
  }
  label 'process_long'
  //when: params.vireo.run_gtmatch_aposteriori

  input:
    tuple val(pool_id), path(vireo_gt_vcf)
    tuple path(ref_gt_vcf), path(ref_gt_csi)

  output:
    path("${donor_assignment_csv}", emit: donor_match_table)
    path("${donor_scores_csv}", emit: donor_score_table)
    path("${gt_check_output_txt}", emit: gtcheck_out)

  script:
    gt_check_output_txt = "${pool_id}_gtcheck.txt"
    score_output_prfx = "${pool_id}"
    donor_assignment_csv = "${score_output_prfx}_gtcheck_donor_assignments.csv"
    donor_scores_csv = "${score_output_prfx}_gtcheck_score_table.csv"
  """
    # fix header of vireo VCF
    #tabix -p vcf ${ref_gt_vcf}
    bcftools view -h ${vireo_gt_vcf} > ${pool_id}_header.txt
    # bcftools view -Ov ${vireo_gt_vcf} > viewed.vcf
    sed -i '/^##fileformat=VCFv.*/a ##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">' ${pool_id}_header.txt
    bcftools reheader -h ${pool_id}_header.txt -o ${pool_id}_GT_donors.vireo.headfix.vcf.gz ${vireo_gt_vcf}

    # sort and index vireo VCF file (bcftools sort bails out with an error)
    bcftools view ${pool_id}_GT_donors.vireo.headfix.vcf.gz | \
      awk '\$1 ~ /^#/ {print \$0;next} {printf"chr%s",\$0 | "sort -k1,1V -k2,2n"}' |
      bcftools view -Oz -o ${pool_id}_chr_GT_donors.vireo.srt.vcf.gz -
    tabix -p vcf ${pool_id}_chr_GT_donors.vireo.srt.vcf.gz

    bcftools gtcheck -g ${ref_gt_vcf} ${pool_id}_chr_GT_donors.vireo.srt.vcf.gz > ${gt_check_output_txt}

    # generate assignment and score tables
    gtcheck_assign.py ${gt_check_output_txt} ${score_output_prfx}

    # generate plots of score density distribution
  """
}
