

def random_hex(n) {
  Long.toUnsignedString(new Random().nextLong(), n).toUpperCase()
}
include {collect_file as collect_file1;
        collect_file as collect_file2;
        collect_file as collect_file3;
        collect_file as collect_file4;
        collect_file as collect_file5;
        collect_file as collect_file6;
        collect_file as collect_file7;
        collect_file as collect_file8} from "$projectDir/modules/local/collect_file/main"


process VACUTAINER_TO_DONOR_ID {
  tag "${study_label}.${pool_id}"
  label 'process_tiny'

  input:
    tuple val(pool_id), val(comma_separated_list_of_vacutainer_ids), val(study_label), path(conversion_file)

  output:
    tuple val(study_label), val(pool_id), path(file_of_donor_ids), emit: study_pool_donorfil optional true

  script:
  file_of_donor_ids = "${study_label}.${pool_id}.donor_ids.txt"
  """
    echo "study_label: ${study_label}"
    echo "pool_id: ${pool_id}"
    echo "list_of_vacutainer_ids: ${comma_separated_list_of_vacutainer_ids}"
    echo "${file_of_donor_ids}"
    vacutainer_to_donor_id.py ${conversion_file} ${comma_separated_list_of_vacutainer_ids} ${file_of_donor_ids}

    # remove file if empty so as to emit no output and stop downstream processes here
    rv=(\$(wc -l ${file_of_donor_ids}))
    if [ \${rv[0]} < 1 ]; then
      rm \${file_of_donor_ids}
    fi
  """
}

process FETCH_DONOR_IDS_FROM_VCF {
  tag "${study_label}.${study_vcf}"

  // return a list of donors from a VCF file
  label 'process_tiny'

  input:
    tuple val(study_label), path(study_vcf)

  output:
    tuple val(study_label), path(vcf_donor_list_file), emit: study_vcf_donor_list
    path "versions.yml", emit: versions

  script:
  vcf_donor_list_file = "${study_label}.vcf_donor_list.txt"
  """
  bcftools query -l ${study_vcf} > ${vcf_donor_list_file}

  cat <<-END_VERSIONS > versions.yml
  "${task.process}":
      bcftools: \$(bcftools --version 2>&1 | head -n1 | sed 's/^.*bcftools //; s/ .*\$//')
  END_VERSIONS
  """
}

process CHECK_DONORS_IN_VCF_HEADER {
  // look up donor ids in VCF header and return a table of <donor_id>,<is_present[Y/N]>
  tag "${pool_id}.${study_label}"
  label 'process_tiny'

  input:
    tuple val(study_label), val(pool_id), path(txt_file_pool_donor_list), path(list_of_vcf_donors_file)

  output:
    tuple val(study_label), val(pool_id), path(tsv_table_of_checked_donor_ids), emit:study_pool_donorfil optional true

  script:
  tsv_table_of_checked_donor_ids = "${pool_id}.${study_label}.donor_ids_checked.tsv"
  """
  echo "study_label: ${study_label}"
  echo "pool_id: ${pool_id}"
  echo "txt_file_pool_donor_list: ${txt_file_pool_donor_list}"
  echo "list_of_vcf_donors_file: ${list_of_vcf_donors_file}"
  echo "output: ${tsv_table_of_checked_donor_ids}"

  find_pooled_donor_ids_in_vcf.py ${list_of_vcf_donors_file} ${txt_file_pool_donor_list} ${tsv_table_of_checked_donor_ids}

  rv=(\$(wc -l ${tsv_table_of_checked_donor_ids}))
  if [ \${rv[0]} < 1 ]; then
    rm \${tsv_table_of_checked_donor_ids}
  fi
  """
}

process SELECT_DONOR_GENOTYPES_FROM_VCF {
  label 'process_tiny'
  if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
      container "${params.yascp_container}"
  } else {
      container "${params.yascp_container_docker}"
  }

  input:
    tuple val(study_label), val(pool_id), path(donor_table), path(study_vcf), path(study_vcf_index)

  output:
    tuple val(study_label), val(pool_id), path(pool_study_bcfgz), emit: study_pool_bcfgz
    path "versions.yml", emit: versions

  script:
  pool_study_bcfgz = "${pool_id}.${study_vcf}.bcf.gz"
  """
    awk 'NR>1 && \$2 !~/^N\$/ {print \$1}' ${donor_table} > donors.lst
    bcftools view -S donors.lst -Ob -o ${pool_study_bcfgz} ${study_vcf}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bcftools: \$(bcftools --version 2>&1 | head -n1 | sed 's/^.*bcftools //; s/ .*\$//')
    END_VERSIONS
  """
}

process CONCAT_STUDY_VCFS {
  label 'process_small'

  if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
      container "${params.yascp_container}"
  } else {
      container "${params.yascp_container_docker}"
  }

  input:
    val(study_label), tuple val(pool_id), path( study_vcf_files )

  output:
    tuple val(study_label), val(pool_id), path(pool_study_bcfgz), emit: study_pool_bcfgz
    path "versions.yml", emit: versions

  script:
  pool_study_bcfgz = "${pool_id}.${study_label}.bcf.gz"
  """
    cat "${study_vcf_files}" > ./fofn_vcfs.txt
    bcftools concat --threads ${task.threads} -f ./fofn_vcfs.txt -Ob -o ${pool_study_bcfgz}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bcftools: \$(bcftools --version 2>&1 | head -n1 | sed 's/^.*bcftools //; s/ .*\$//')
    END_VERSIONS
  """
}

process SUBSET_GENOTYPE {
    tag "${samplename}.${sample_subset_file}"
    label 'process_medium'
    publishDir "${params.outdir}/preprocessing/subset_genotypes/", mode: "${params.copy_mode}", pattern: "${samplename}.${sample_subset_file}.subset.vcf.gz"


    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "${params.yascp_container}"
    } else {
        container "${params.yascp_container_docker}"
    }

    input:
    tuple val(samplename), path(donor_vcf),path(donor_vcf_csi), val(sample_subset_file)


    output:
    tuple val(samplename), path("${samplename}.subset.vcf.gz"),path("${samplename}.subset.vcf.gz.csi"), emit: samplename_subsetvcf
    path "versions.yml", emit: versions
    
    script:
    """
        echo ${sample_subset_file}
        #tabix -p vcf ${donor_vcf} || echo 'not typical VCF'
        bcftools view ${donor_vcf} -s ${sample_subset_file} -Oz -o ${samplename}.subset.vcf.gz
        bcftools index ${samplename}.subset.vcf.gz
        rm ${donor_vcf}.tbi || echo 'not typical VCF'

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            bcftools: \$(bcftools --version 2>&1 | head -n1 | sed 's/^.*bcftools //; s/ .*\$//')
        END_VERSIONS
    """
}


process SUBSET_GENOTYPE2 {
    tag "${cohort}___${samplename}"
    label 'process_low'

    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "${params.yascp_container}"
    } else {
        container "${params.yascp_container_docker}"
    }

    input:
      tuple val(samplename), val(sample_subset_file),val(cohort),path(donor_vcf),path(donor_vcf_csi)


    output:
      tuple val("${cohort}___${samplename}"), path("${donor_vcf}_subset.vcf.gz"),path("${donor_vcf}_subset.vcf.gz.csi"), emit: subset_vcf_file optional true
      path('*_mapping.tsv'), emit: mapping optional true
      path "versions.yml", emit: versions
    script:
      if (params.genotype_phenotype_mapping_file!=''){
        g_p_map = " -b ${params.genotype_phenotype_mapping_file}"
      }else{
        g_p_map = " "
      }

    """
        echo "${cohort}___${samplename}"
        echo *_${donor_vcf}_subset.vcf.gz
        bcftools query -l ${donor_vcf} > samples.tsv
        extract_overlaps.py -vs samples.tsv ${g_p_map} -s ${sample_subset_file} -o sample_file.tsv || touch sample_file.tsv && echo 'no input, as a result of samples missing from the file'
        if [ \$(cat sample_file.tsv | wc -l) -gt 0 ]; then
          bcftools view ${donor_vcf} -S sample_file.tsv -Oz -o ${donor_vcf}_subset.vcf.gz 
          bcftools index ${donor_vcf}_subset.vcf.gz
        else
          echo 'no'
        fi
        rm samples.tsv sample_file.tsv

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            bcftools: \$(bcftools --version 2>&1 | head -n1 | sed 's/^.*bcftools //; s/ .*\$//')
        END_VERSIONS
    """
}


process JOIN_CHROMOSOMES{
    tag "${samplename}"
    label 'process_medium'
    publishDir "${params.outdir}/preprocessing/subset_genotypes/", mode: "${params.copy_mode}", pattern: "${samplename}.${sample_subset_file}.subset.vcf.gz"

    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "${params.yascp_container}"
    } else {
        container "${params.yascp_container_docker}"
    }

    input:
      tuple val(samplename), path(study_vcf_files),path(study_vcf_csi_files)
      path(genome)


    output:
      tuple val(s2), path("*_out.bcf.gz"),path("*_out.bcf.gz.csi"), emit: joined_chromosomes_per_studytrance
      path "versions.yml", emit: versions

    script:
      s1 = samplename.split('___')[0]
      s2 = samplename.split('___')[1]

      """
        vcf_name=\$(python ${projectDir}/bin/random_id.py)
        fofn_input_subset.sh "${study_vcf_files}"
        bcftools concat --threads ${task.threads} -f ./fofn_vcfs.txt -Ob -o pre_\${vcf_name}.bcf.gz

        bcftools view pre_\${vcf_name}.bcf.gz | awk '{gsub(/^chr/,""); print}' | awk '{gsub(/ID=chr/,"ID="); print}' > no_prefix_pre_\${vcf_name}.vcf
        bgzip no_prefix_pre_\${vcf_name}.vcf
        bcftools index no_prefix_pre_\${vcf_name}.vcf.gz
        bcftools view ${params.bcf_viewfilters} no_prefix_pre_\${vcf_name}.vcf.gz -Oz -o no_prefix2_pre_\${vcf_name}.vcf.gz
        bcftools index no_prefix2_pre_\${vcf_name}.vcf.gz
        #bcftools index pre_\${vcf_name}.bcf.gz
        bcftools +fixref no_prefix2_pre_\${vcf_name}.vcf.gz -Ob -o fix_ref_\${vcf_name}_out.bcf.gz -- -d -f ${genome}/genome.fa -m flip-all
        #ln -s pre_\${vcf_name}.bcf.gz \${vcf_name}_out.bcf.gz
        #bcftools +fill-tags fix_ref_\${vcf_name}_out.bcf.gz -Ob -o \${vcf_name}_out.bcf.gz
        bcftools annotate -x INFO fix_ref_\${vcf_name}_out.bcf.gz -Ob -o \${vcf_name}_out.bcf.gz
        rm -r pre_* fix_ref_* no_prefix_*

        bcftools index \${vcf_name}_out.bcf.gz

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            bcftools: \$(bcftools --version 2>&1 | head -n1 | sed 's/^.*bcftools //; s/ .*\$//')
            bgzip: \$(echo \$(bgzip -h 2>&1) | head -n 1 | sed 's/^.*(htslib) //; s/ .*\$//')
        END_VERSIONS
      """
}

process RESOLVE_POOL_VCFS{
    tag "${samplename}"
    publishDir "${params.outdir}/preprocessing/subset_genotypes", mode: "${params.copy_mode}",
    saveAs: {filename ->
          if (filename.contains("AllExpectedGT")) {
            if (filename.contains("Data_Pipeline___")) {
              null
            } else if (filename.contains("Data_User___"))  {
                null
            }else{
              filename
            }
          } else {
            null
          }
        },
      overwrite: "true"

    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "${params.yascp_container}"
    } else {
        container "${params.yascp_container_docker}"
    }
    input:
      tuple val(samplename), path(vcf),path(vcf_csi)
      val(mode)
    output:
      path('Data_Pipeline___*'), emit: pipeline_data
      path('Genotype___*'), emit: genotype_folder
      path('Data_User___*'), emit: user_data
      // If more than one pool is using the same genotype it is pointless to emit it many times. Hence we produce a vcf pointer files which indicate which pool uses which genotype.
    script:
      vcf=vcf[0]
      """
        pool_panel.py --mode ${mode} --vcf ${vcf} --pool_ids ${samplename} --mode ${mode}
        
      """
}

process JOIN_STUDIES_MERGE{
    tag "${samplename}"
    label 'process_medium_memory'

    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "${params.yascp_container}"
    } else {
        container "${params.yascp_container_docker}"
    }

    input:
      tuple val(samplename), path(study_vcf_files),path(study_vcf_csi_files)
      val(mode)
      val(mode2)


    output:
      tuple val(samplename), path("*_out.vcf.gz"),path("*_out.vcf.gz.csi"), emit: merged_expected_genotypes
      path("*_out.vcf.gz",emit:study_merged_vcf)
      path "versions.yml", emit: versions
    script:

        
        
        if (params.just_overlapping_positions_for_study_merge){
          cmd__run = "overlapping_positions_vcfs.py -vcfs '${study_vcf_files}'"
          cmd="bcftools view -R Bed_File_record.bed pre_${mode}_${mode2}_\${vcf_name}__vcf.vcf.gz -Oz -o ${mode}_${mode2}_\${vcf_name}_out.vcf.gz && bcftools index ${mode}_${mode2}_\${vcf_name}_out.vcf.gz"
        }else{
          cmd__run = ''
          cmd="mv pre_${mode}_${mode2}_\${vcf_name}__vcf.vcf.gz ${mode}_${mode2}_\${vcf_name}_out.vcf.gz && mv pre_${mode}_${mode2}_\${vcf_name}__vcf.vcf.gz.csi ${mode}_${mode2}_\${vcf_name}_out.vcf.gz.csi"
        }

      """
        vcf_name=\$(python ${projectDir}/bin/random_id.py)
       
        
        fofn_input_subset.sh "${study_vcf_files}"
        if [ \$(cat fofn_vcfs.txt | wc -l) -gt 1 ]; then
            echo 'yes'
            ${cmd__run}
            bcftools merge -file-list ${study_vcf_files} -Ou | bcftools sort -T \$PWD -Oz -o pre_${mode}_${mode2}_\${vcf_name}__vcf.vcf.gz
            bcftools index pre_${mode}_${mode2}_\${vcf_name}__vcf.vcf.gz
            ${cmd}
        else
          echo 'no'
          bcftools view ${study_vcf_files} | bcftools sort -T \$PWD -Oz -o ${mode}_${mode2}_\${vcf_name}_out.vcf.gz
          bcftools index ${mode}_${mode2}_\${vcf_name}_out.vcf.gz 
          
        fi
     
        rm -r pre_* || echo 'nothing to remove'

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            bcftools: \$(bcftools --version 2>&1 | head -n1 | sed 's/^.*bcftools //; s/ .*\$//')
        END_VERSIONS
      """
}


workflow SUBSET_WORKF{
  take:
    ch_ref_vcf
    donors_in_pools
    mode
    genome
  main:
      donors_in_pools.combine(ch_ref_vcf).unique().set{all_GT_pannels_and_pools}
      all_GT_pannels_and_pools.map { row -> tuple("${row[1]}:${row[2]}:${row[3]}:${row[4]}",row[0],row[1],row[2],row[3],row[4]) }.set { combined_pool_subset }

      combined_pool_subset.groupTuple(by: 0).set{grouped_chrs_poolComps}

      if (!params.genotype_input.subset_vireo_genotypes && mode=='AllExpectedGT'){
        // Here we have joned the shards of the same cohort together without subsetting down to the individuals since we want to use all the genotypes as an input file in the vcf.
        // While this will work, if the same cohort is used multiple times its better to provide already merged genotype, to avoid all the shard processing every time.
        grouped_chrs_poolComps.map { row -> tuple( row[0], "${row[0]}".split(':')[0],  "${row[0]}".split(':')[1],  "${row[0]}".split(':')[2],  "${row[0]}".split(':')[3],row[1].join('::') ) }.set { combined_pool_subset } // Here we have all the pools that contain the same donor compositions associated with each of the shards
        combined_pool_subset.map{row -> tuple( "${row[2]}____${row[5]}", row[3],row[4])}.set{chromosome_vcfs_per_studypool1}
        chromosome_vcfs_per_studypool1.groupTuple().set{chromosome_vcfs_per_studypool}
      }else{
        // Here we subset down each of the unique donor pools in each of the shards.
        // We utilise this mode to also check each of the shards against expected and gt matched genotypes. 
        // In case where we do not have any info re what is expected but we have run it in genotype avare mode, we do not perform the IBD calculations.
        grouped_chrs_poolComps.map { row -> tuple( row[1].join('::'), "${row[0]}".split(':')[0],  "${row[0]}".split(':')[1],  "${row[0]}".split(':')[2],  "${row[0]}".split(':')[3] ) }.set { combined_pool_subset } // Here we have all the pools that contain the same donor compositions associated with each of the shards
        grouped_chrs_poolComps.map { row -> tuple( row[0], row[1]) }.set { pools_utilising_same_subset } // Here we have a mapping file of which pools should use which genotypes.
        SUBSET_GENOTYPE2(combined_pool_subset)

        SUBSET_GENOTYPE2.out.subset_vcf_file.unique().groupTuple().set{chromosome_vcfs_per_studypool}

      }

      // Combining all the chromosomes per pool
      // Now we combine all the chromosomes together for all the unique pool compositions.
      JOIN_CHROMOSOMES(chromosome_vcfs_per_studypool,genome)
      JOIN_CHROMOSOMES.out.joined_chromosomes_per_studytrance.unique().groupTuple().set{study_vcfs_per_pool}
      // study_vcfs_per_pool.subscribe {println "study_vcfs_per_pool:= ${it}\n"}

      // Since user are capable in providing multiple different study vcfs these need to be merged together per unique pool composition to perform the internal IBD checks.
      JOIN_STUDIES_MERGE(study_vcfs_per_pool,'Study_Merge',mode)
      JOIN_STUDIES_MERGE.out.merged_expected_genotypes.set{merged_expected_genotypes}

      // After merging studies per unique pool compositions we resolve the matches back to the each of the Pools so that the IBD and Vireo can use the correct genotypes as the inputs and publish these in the correct folder.
      RESOLVE_POOL_VCFS(JOIN_STUDIES_MERGE.out.merged_expected_genotypes,mode)
      pools_panels = RESOLVE_POOL_VCFS.out.pipeline_data

      if (mode=='AllExpectedGT'){
        collect_file1(RESOLVE_POOL_VCFS.out.user_data.collect(),"Genotypes_all_pools.tsv",params.outdir+'/preprocessing/subset_genotypes',1,'')
      }
      pools_panels.splitCsv(header: true, sep: '\t').map { row -> tuple(row['Pool_id'], file(row.vcf), file(row.vcf_csi)) }
                .set{merged_expected_genotypes}
      pools_panels.splitCsv(header: true, sep: '\t').map { row -> tuple(row['Pool_id'], file(row.vcf)) }
                .set{samplename_subsetvcf_ibd}
      // If we are not subsetting anything down to the expected donors and are planning to use these genotypes multiple times we should just provide a single cohort genotype file as an input.
      study_merged_vcf = JOIN_STUDIES_MERGE.out.study_merged_vcf
      

  emit:
    merged_expected_genotypes
    study_merged_vcf
    samplename_subsetvcf_ibd
}
