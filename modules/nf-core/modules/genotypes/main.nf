process MERGE_GENOTYPES_IN_ONE_VCF_IDX_PAN{

    label 'process_medium'
    publishDir  path: "${params.outdir}/${mode}_genotypes/${pn1}",
          mode: "${params.copy_mode}",
          overwrite: "true"
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        // println "container: /software/hgi/containers/wtsihgi-nf_genotype_match-1.0.sif\n"
        container "${params.nf_yascp_htstools_container}"
    } else {
        container "mercury/wtsihgi-nf_yascp_htstools-1.1"
    }

    input:
       tuple val(panel), path(vireo_gt_vcf), path(vireo_gt_vcf_csi),path(barcodes)
       val(mode)

    output:
       tuple  val(pn1), path("${mode}.${panel}.vcf.gz"),path("${mode}.${panel}.vcf.gz.csi"),path(barcodes), emit: gt_pool
      //  path("${mode}.${panel}.vcf.gz"), emit: study_merged_vcf optional true
      // here we want to make it look like its a vireo output file
    script:
        def pan = "${panel}".tokenize('.')
        pn1 = pan[0]
    """
      echo ${pn1}
      fofn_input_subset.sh "${vireo_gt_vcf}"

      if [ \$(cat fofn_vcfs.txt | wc -l) -gt 1 ]; then
          echo 'yes'
          bcftools concat -f fofn_vcfs.txt -Ou | bcftools sort -T \$PWD -Oz -o ${mode}.${panel}.vcf.gz
          bcftools index ${mode}.${panel}.vcf.gz
      else
        echo 'no'
        bcftools view ${vireo_gt_vcf} | bcftools sort -T \$PWD -Oz -o ${mode}.${panel}.vcf.gz
        bcftools index ${mode}.${panel}.vcf.gz
      fi
    """

}


process MERGE_GENOTYPES_IN_ONE_VCF_FREEBAYES{

    label 'process_medium'
    publishDir  path: "${params.outdir}/${mode}_genotypes",
              saveAs: {filename ->
                    if (filename.endsWith("vireo_${panel}")) {
                        null
                    } else{
                        filename
                    }
                },
          mode: "${params.copy_mode}",
          overwrite: "true"

    publishDir  path: "${params.outdir}/deconvolution/vireo",
          saveAs: {filename ->
                    if (filename.endsWith("vireo_${panel}")) {
                        filename
                    } else{
                        null
                    }
                },
          mode: "${params.copy_mode}",
          overwrite: "true"
    
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        // println "container: /software/hgi/containers/wtsihgi-nf_genotype_match-1.0.sif\n"
        container "${params.nf_yascp_htstools_container}"
    } else {
        container "mercury/wtsihgi-nf_yascp_htstools-1.1"
    }

    input:
       tuple val(panel), path(vireo_gt_vcf), path(vireo_gt_vcf_csi),path(barcodes)
       val(mode)

    output:
      tuple  val(panel), path("${mode}.${panel}.vcf.gz"),path("${mode}.${panel}.vcf.gz.csi"), emit: gt_pool
      //  path("${mode}.${panel}.vcf.gz"), emit: study_merged_vcf optional true
      path("vireo_${panel}"), emit: vir_input
    script:

    """
      echo ${panel}
      fofn_input_subset.sh "${vireo_gt_vcf}"

      if [ \$(cat fofn_vcfs.txt | wc -l) -gt 1 ]; then
          echo 'yes'
          bcftools merge --force-samples -file-list ${vireo_gt_vcf}  -Ou | bcftools sort -T \$PWD -Oz -o ${mode}.${panel}.vcf.gz
          bcftools index ${mode}.${panel}.vcf.gz
      else
        echo 'no'
        bcftools view ${vireo_gt_vcf} | bcftools sort -T \$PWD -Oz -o ${mode}.${panel}.vcf.gz
        bcftools index ${mode}.${panel}.vcf.gz
      fi

      process_barcodes.sh
      mkdir vireo_${panel}
      cd vireo_${panel} && ln -s ../${mode}.${panel}.vcf.gz ./GT_donors.vireo.vcf.gz && mv ../donor_ids.tsv ./
    """

}

process MERGE_GENOTYPES_IN_ONE_VCF{

    label 'process_medium'
    publishDir  path: "${params.outdir}/${mode}_genotypes/",
          mode: "${params.copy_mode}",
          overwrite: "true"
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        // println "container: /software/hgi/containers/wtsihgi-nf_genotype_match-1.0.sif\n"
        container "${params.nf_yascp_htstools_container}"
    } else {
        container "mercury/wtsihgi-nf_yascp_htstools-1.1"
    }

    input:
       path(vireo_gt_vcf)
       val(mode)

    output:
       tuple path("${mode}_merged_vcf_file_all_pools.vcf.gz"),path("${mode}_merged_vcf_file_all_pools.vcf.gz.csi"), emit: merged_infered_genotypes
      //  path("${mode}_merged_vcf_file_all_pools.vcf.gz"), emit: study_merged_vcf optional true

    script:

    """
      fofn_input_subset.sh "${vireo_gt_vcf}"

      for VARIABLE in ${vireo_gt_vcf}
      do
          bcftools index -f \$VARIABLE
      done

      if [ \$(cat fofn_vcfs.txt | wc -l) -gt 1 ]; then
          echo 'yes'
          bcftools merge --force-samples -i MAF:join -file-list ${vireo_gt_vcf} -Ou | bcftools sort -T \$PWD -Oz -o ${mode}_merged_vcf_file_all_pools.vcf.gz
          bcftools index ${mode}_merged_vcf_file_all_pools.vcf.gz
      else
        echo 'no'
        bcftools view ${vireo_gt_vcf} | bcftools sort -T \$PWD -Oz -o ${mode}_merged_vcf_file_all_pools.vcf.gz
        bcftools index ${mode}_merged_vcf_file_all_pools.vcf.gz
        
      fi
    """

}


process VIREO_ADD_SAMPLE_PREFIX{

    tag "${pool_id}"
    label 'process_low'

    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        // println "container: /software/hgi/containers/wtsihgi-nf_genotype_match-1.0.sif\n"
        container "${params.nf_yascp_htstools_container}"
    } else {
        container "mercury/wtsihgi-nf_yascp_htstools-1.1"
    }

    input:
      tuple val(pool_id), path(vireo_gt_vcf)

    output:
      path("prefix_${vireo_fixed_vcf}"), emit: infered_vcf

    script:
      sorted_vcf = "${pool_id}_vireo_srt.vcf.gz"
      vireo_fixed_vcf = "${pool_id}_pool_headfix_vireo.vcf.gz"
    """
      bcftools query -l ${vireo_gt_vcf} | awk '\$0=""\$0" ${pool_id}_"\$0' > replacement_assignments.tsv
      bcftools reheader --samples replacement_assignments.tsv -o prefix_${vireo_fixed_vcf} ${vireo_gt_vcf}
    """
}

process VIREO_GT_FIX_HEADER
{
  tag "${pool_id}"
  publishDir  path: "${params.outdir}/infered_genotypes/${pool_id}/",
        mode: "${params.copy_mode}",
        overwrite: "true"
  if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
      // println "container: /software/hgi/containers/wtsihgi-nf_genotype_match-1.0.sif\n"
      container "${params.nf_yascp_celltypist}"
  } else {
      container "mercury/wtsihgi-nf_yascp_htstools-1.1"
  }

  label 'process_low'

  input:
    tuple val(pool_id), path(vireo_gt_vcf)
    path(genome)

  output:
    tuple val(pool_id), path("${vireo_fixed_vcf}"), path("${vireo_fixed_vcf}.tbi"), emit: gt_pool
    tuple val(pool_id), path("${vireo_fixed_vcf}"), emit: infered_vcf

  script:
  sorted_vcf = "${pool_id}_vireo_srt.vcf.gz"
  vireo_fixed_vcf = "${pool_id}_headfix_vireo.vcf.gz"


  """
    # fix header of vireo VCF
    bcftools view -h ${vireo_gt_vcf} > init_head.txt
    sed -i '/^##fileformat=VCFv.*/a ##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">' init_head.txt
    head -n -1 init_head.txt > header.txt
    echo '##INFO=<ID=AD,Number=A,Type=Integer,Description="alternative allele  (variant-by-cell) of reads">' >> header.txt
    echo '##INFO=<ID=DP,Number=1,Type=Integer,Description="depth UMIs for each variant in each cell">' >> header.txt
    echo '##INFO=<ID=PL,Number=1,Type=Integer,Description="depth UMIs for each variant in each cell">' >> header.txt
    echo '##INFO=<ID=OTH,Number=1,Type=Integer,Description="????">' >> header.txt
    echo '##FORMAT=<ID=PL,Number=G,Type=Integer,Description="???">' >> header.txt
    echo '##FORMAT=<ID=AD,Number=G,Type=Integer,Description="????n">' >> header.txt
    echo '##FORMAT=<ID=DP,Number=G,Type=Integer,Description="????n">' >> header.txt
    tail -n1 init_head.txt >> header.txt

    # sort VCF file (bcftools sort bails out with an error)
    bcftools view ${vireo_gt_vcf} | \
    awk '\$1 ~ /^#/ {print \$0;next} {print \$0 | "sort -k1,1V -k2,2n"}' | \
    bcftools view -Oz -o ${sorted_vcf} -

    bcftools reheader -h header.txt ${sorted_vcf} | \
    bcftools view -Oz -o pre_${vireo_fixed_vcf}
    tabix -p vcf pre_${vireo_fixed_vcf}
    bcftools +fixref pre_${vireo_fixed_vcf} -Oz -o ${vireo_fixed_vcf} -- -d -f ${genome}/genome.fa -m flip-all
    tabix -p vcf ${vireo_fixed_vcf}


  """
}
process REPLACE_GT_DONOR_ID2{
    tag "${samplename}"
    publishDir  path: "${params.outdir}/deconvolution/vireo_gt_fix/${samplename}/",
          pattern: "GT_replace_*",
          mode: "${params.copy_mode}",
          overwrite: "true"
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "${params.scrna_deconvolution}"
    } else {
        container "mercury/scrna_deconvolution:62bd56a"
    }

  label 'process_low'

  input:
    tuple val(samplename), path(gt_donors), path(vireo_sample_summary),path(vireo___exp_sample_summary),path(vireo__donor_ids),path(vcf_file),path(donor_gt_csi),val(mode)

  output:
    tuple val(samplename), path("GT_replace_donor_ids_${mode}.tsv"), emit: sample_donor_ids
    tuple val(samplename), path("GT_replace_GT_donors.vireo_${mode}.vcf.gz"), path(vcf_file),path(donor_gt_csi), emit: sample_donor_vcf
    tuple val(samplename), path("GT_replace_GT_donors.vireo_${mode}.vcf.gz"), emit: infered_vcf
    path("GT_replace_${samplename}.sample_summary_${mode}.txt"), emit: sample_summary_tsv
    path("GT_replace_${samplename}__exp.sample_summary_${mode}.txt"), emit: sample__exp_summary_tsv
    path("GT_replace_${samplename}_assignments_${mode}.tsv"), emit: assignments
    tuple  val(samplename), path("GT_replace_GT_donors.vireo_${mode}.vcf.gz"), path("GT_replace_${samplename}.sample_summary_${mode}.txt"),path("GT_replace_${samplename}__exp.sample_summary_${mode}.txt"),path("GT_replace_donor_ids_${mode}.tsv"),path(vcf_file),path(donor_gt_csi), emit: all_required_data
    tuple val(samplename), path("GT_replace_donor_ids_${mode}.tsv"), emit: cell_assignments
    path("vireo_gt_rep_${samplename}"), emit: output_dir
  script:

    in=""

    """
      bcftools query -l ${gt_donors} > ${mode}_donors_in_vcf.tsv
      replace_donors.py -id ${samplename} ${in} --input_file "${params.input_data_table}" -m ${mode}
      bcftools view ${gt_donors} | bcftools reheader --samples replacement_assignments_${mode}.tsv -o GT_replace_GT_donors.vireo_${mode}.vcf.gz
      mkdir -p "vireo_gt_rep_${samplename}"
      ln -sf "\$(realpath "GT_replace_GT_donors.vireo_${mode}.vcf.gz")" "vireo_gt_rep_${samplename}/GT_replace_GT_donors.vireo_${mode}.vcf.gz"
    """
}

process ENHANCE_STATS_GT_MATCH{

    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "${params.scrna_deconvolution}"
    } else {
        container "mercury/scrna_deconvolution:62bd56a"
    }
  tag "${samplename}"
  publishDir  path: "${params.outdir}/gtmatch/${samplename}",
          mode: "${params.copy_mode}",
          overwrite: "true"

  label 'process_medium'

  input:
    tuple val(samplename), path(enhancement_file)
    path(input_data_table)
  output:

    path("GT_replace_${enhancement_file}"), emit: assignments
    
  script:
    if(params.genotype_phenotype_mapping_file==''){
      in=""
    }else if (params.use_phenotype_ids_for_gt_match){
      in="--genotype_phenotype_mapping ${params.genotype_phenotype_mapping_file}"
      // in=""
    }else{
      in=""
    }

    // if(params.input_data_table==''){
    //   in_f = ""
    // }else{
    //   in_f = 
    // }

    """
      enhance_stats.py -id ${samplename} -dm ${enhancement_file} ${in} --input_file '${input_data_table}' -m ${params.genotype_input.vireo_with_gt}
    """
}

process GT_MATCH_POOL_IBD
{
  tag "${pool_id}_ibd"
  label 'process_small'
  publishDir  path: "${params.outdir}/gtmatch/${pool_id}",
          mode: "${params.copy_mode}",
          overwrite: "true"

  if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
      container "${params.nf_scrna_qc_v3_container}"
  } else {
      container "mercury/nf_qc_scrna:v2"
  }

  

  input:
    tuple val(pool_id), path(vireo_gt_vcf)
    val(mode)
    val(mode2)

  output:
    tuple val(pool_id),path("*_${pool_id}.genome*"), emit:plink_ibd optional true

  script:
    """
      #bcftools +prune -m 0.2 -w 50 ${vireo_gt_vcf} -Ov -o pruned_${vireo_gt_vcf}
      #plink --vcf ${vireo_gt_vcf} --indep-pairwise 50 5 0.2 --out all2 --make-bed --double-id
      #plink --bfile all2 --extract all2.prune.in --out pruned --export vcf
      plink --vcf ${vireo_gt_vcf} --genome unbounded --const-fid dummy --out ${mode2}_${mode}_${pool_id} || echo 'single individual pool, cant calculate IBD'
      #rm all*
    """
}

process GT_MATCH_POOL_AGAINST_PANEL
{
  tag "${pool_id}_vs_${panel_id}"

  if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
      // println "container: /software/hgi/containers/wtsihgi-nf_genotype_match-1.0.sif\n"
      container "${params.nf_yascp_htstools_container}"
  } else {
      container "mercury/wtsihgi-nf_yascp_htstools-1.1"
  }

  label 'process_tiny'
  //when: params.vireo.run_gtmatch_aposteriori

  input:
    tuple val(pool_id), path(vireo_gt_vcf), path(vireo_gt_tbi), val(panel_id), path(ref_gt_vcf), path(ref_gt_csi)

  output:
    tuple val(pool_panel_id), path("${gt_check_output_txt}"), emit:gtcheck_results

  script:
  pool_panel_id = "pool_${pool_id}_panel_${panel_id}"
  panel_filnam = "${ref_gt_vcf}" - (~/\.[bv]cf(\.gz)?$/)
  gt_check_output_txt = "${pool_id}_gtcheck_${panel_filnam}.txt"
  """
    bcftools gtcheck --no-HWE-prob -g ${ref_gt_vcf} ${vireo_gt_vcf} > ${gt_check_output_txt}
  """
}


process PREPROCESS_GENOTYPES
{
  tag "${pool_id}_vs_${panel_id}"

  if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
      // println "container: /software/hgi/containers/wtsihgi-nf_genotype_match-1.0.sif\n"
      container "https://yascp.cog.sanger.ac.uk/public/singularity_images/wtsihgi-nf_yascp_htstools-1.1.sif"
  } else {
      container "mercury/wtsihgi-nf_yascp_htstools-1.1"
  }

  label 'process_tiny'
  input:
    tuple val(pool_id), path(ref_gt_vcf), path(ref_gt_csi)

  output:
    tuple val(pool_id), path("renamed_*.vcf.gz"), path("renamed_*.vcf.gz.csi")

  script:

  """
    renamed_vcf_basename=\$(basename "${ref_gt_vcf}" | sed -E 's/\\.(vcf|bcf)(\\.gz)?\$//')
    renamed_vcf="renamed_\${renamed_vcf_basename}.vcf" 

    # Check if the VCF file has chromosome prefixes
    STR=\$(bcftools index -s ${ref_gt_vcf} | cut -f1 | head -n1 || echo "no_chr")
    SUB='chr'
    if [[ "\$STR" == *"\$SUB"* ]]; then
      # Remove 'chr' prefix and re-save with the 'renamed_' prefix
      zcat "${ref_gt_vcf}" | awk '{gsub(/^chr/,""); print}' | awk '{gsub(/ID=chr/,"ID="); print}' > "\${renamed_vcf}"
      bgzip "\${renamed_vcf}"  # bgzip will add .gz automatically
      bcftools index "\${renamed_vcf}.gz"
    else
      # Create symbolic links with 'renamed_' prefix
      ln -s "${ref_gt_vcf}" "renamed_\${renamed_vcf_basename}.vcf.gz"
      ln -s "${ref_gt_csi}" "renamed_\${renamed_vcf_basename}.vcf.gz.csi"
    fi
  """
}

process ASSIGN_DONOR_FROM_PANEL
{
  // sum gtcheck discrepancy scores from multiple ouputput files of the same panel
  tag "${pool_panel_id}"
  label 'process_medium'
  publishDir  path: "${params.outdir}/gtmatch/${pool_id}",
          pattern: "*.csv",
          mode: "${params.copy_mode}",
          overwrite: "true"
  if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
      // println "container: /software/hgi/containers/wtsihgi-nf_genotype_match-1.0.sif\n"
      container "${params.nf_scrna_qc_v3_container}"
  } else {
      container "mercury/wtsihgi-nf_genotype_match-1.0"
  }

  input:
    tuple val(pool_panel_id), path(gtcheck_output_files)

  output:
    tuple val(pool_id), path("${assignment_table_out}"), emit: gtcheck_assignments
    path("${score_table_out}", emit: gtcheck_scores)

  

  script:
  (_, pool_id) = ("${pool_panel_id}" =~ /^pool_(\S+)_panel_/)[0]
  score_table_out = "${pool_panel_id}_gtcheck_score_table.csv"
  assignment_table_out = "${pool_panel_id}_gtcheck_donor_assignments.csv"

  """
    gtcheck_assign.py ${pool_panel_id} ${gtcheck_output_files}
  """
}

process ASSIGN_DONOR_OVERALL
{
  // decide final donor assignment across different panels from per-panel donor assignments
  tag "${pool_id}"

  publishDir  path: "${params.outdir}/gtmatch/${pool_id}",
          pattern: "*.csv",
          mode: "${params.copy_mode}",
          overwrite: "true"

  if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
      // println "container: /software/hgi/containers/wtsihgi-nf_genotype_match-1.0.sif\n"
      container "${params.scrna_deconvolution}"
  } else {
      container "mercury/wtsihgi-nf_genotype_match-1.0"
  }

  input:
    tuple val(pool_id), path(gtcheck_assign_files)

  output:
    tuple val(pool_id), path("${donor_assignment_file}"), emit: donor_assignments
    path(stats_assignment_table_out), emit: donor_match_table
    tuple val(pool_id),path(stats_assignment_table_out), emit: donor_match_table_with_pool_id
    path("*.csv")

  label 'process_tiny'

  script:
  donor_assignment_file = "${pool_id}_gt_donor_assignments.csv"
  stats_assignment_table_out = "stats_${pool_id}_gt_donor_assignments.csv"
  """
    gtcheck_assign_summary.py ${donor_assignment_file} ${params.genotype_input.ZSCORE_THRESH} ${params.genotype_input.ZSCORE_DIST_THRESH} ${gtcheck_assign_files}
  """
}




process REPLACE_GT_ASSIGNMENTS_WITH_PHENOTYPE{
  label 'process_low'
  publishDir  path: "${params.outdir}/gtmatch/",
          pattern: "*_assignments.csv",
          mode: "${params.copy_mode}",
          overwrite: "true"

  if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
      container "${params.scrna_deconvolution}"
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

process ENHANCE_STATS_FILE{

  tag "${pool_id}"

  publishDir  path: "${params.outdir}/gtmatch/${pool_id}",
        mode: "${params.copy_mode}",
        overwrite: "true"


  if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
      container "${params.scrna_deconvolution}"
  } else {
      container "mercury/scrna_deconvolution:62bd56a"
  }

  label 'process_small'


  input:
    tuple val(pool_id),path(ibd_table),path(stats_table),val(expected_ids),path(withinn_pool_ibd)
    val(condition)


  output:
    
    tuple val(pool_id),path("PiHAT_Stats_File_${pool_id}.csv"), emit: stats_table_PiHat_enhanced
    path ('Max_PiHAT_For_Expected*'), emit: max_PiHAT_For_Expected optional true
    path('Done.tmp'), emit: done_validation
  script:
    if (params.extra_sample_metadata==''){
      md_inp = ""
    }else{
      md_inp = "-md ${params.extra_sample_metadata}"
    }

    if (params.genotype_phenotype_mapping_file==''){
      mapping=""
    }else{
      mapping="-m ${params.genotype_phenotype_mapping_file} "
    }

    """
      add_PiHat_to_GT_match.py -mt ${stats_table} -ph ${ibd_table} ${mapping} -c ${condition} -e ${expected_ids} -id ${pool_id} -wpi ${withinn_pool_ibd} ${md_inp} || echo 'we dont have expected samples in this cohort'
      echo 'Done' > Done.tmp
    """

}

process ENHANCE_VIREO_METADATA_WITH_DONOR{
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "${params.scrna_deconvolution}"
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



process COMBINE_MATCHES_IN_EXPECTED_FORMAT{
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "${params.scrna_deconvolution}"
    } else {
        container "mercury/scrna_deconvolution:62bd56a"
    }
  label 'process_small'

  input:
    path(stats_files)

  output:
    path('All_Infered_Expected.csv'), emit: all_Infered_Expected

  script:
    if (params.cohorts_to_drop_from_GT_Relatednes_check==''){
      md_in = ""
    }else{
      md_in = "-dr ${params.cohorts_to_drop_from_GT_Relatednes_check}"
    }
    """
      combine_all_GTmatched_in_expected_format.py --files "${stats_files}" ${md_in}
    """
}


workflow MATCH_GT_VIREO {
  take:
    gt_math_pool_against_panel_input


  main:
    // now match genotypes against a panels
    GT_MATCH_POOL_AGAINST_PANEL(gt_math_pool_against_panel_input)

    // group by panel id
    GT_MATCH_POOL_AGAINST_PANEL.out.gtcheck_results.unique()
      .groupTuple()
      .set { gt_check_by_panel }
    

    ASSIGN_DONOR_FROM_PANEL(gt_check_by_panel)
    ASSIGN_DONOR_FROM_PANEL.out.gtcheck_assignments.unique()
      .groupTuple()
      .set{ ch_donor_assign_panel }

    ASSIGN_DONOR_OVERALL(ch_donor_assign_panel)

  emit:
    pool_id_donor_assignments_csv = ASSIGN_DONOR_OVERALL.out.donor_assignments
    donor_match_table = ASSIGN_DONOR_OVERALL.out.donor_match_table
    donor_match_table_with_pool_id = ASSIGN_DONOR_OVERALL.out.donor_match_table_with_pool_id
}
