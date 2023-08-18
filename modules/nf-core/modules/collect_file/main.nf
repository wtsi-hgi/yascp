process collect_file{
  label 'process_tiny'
  if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
    container "https://yascp.cog.sanger.ac.uk/public/singularity_images/wtsihgi_nf_scrna_qc_6bb6af5-2021-12-23-3270149cf265.sif"
    //// container "/lustre/scratch123/hgi/projects/ukbb_scrna/pipelines/singularity_images/wtsihgi_nf_cellbender_v1.2.img"
  } else {
    container "wtsihgi/nf_scrna_qc:6bb6af5"
  }
    // In nf there is a function collectFile - however if you provide a symlinked file directory tusing nf function will overwrite the source instead of replacing the file
    // This snipped is a replication of the function but as a nf module and hence the problem is avoided.
  publishDir "${outpath}/",  mode: "${params.copy_mode}", overwrite: true
  input:
    path(files_to_concentrate)
    val(name)
    val(outpath)
    val(header_size)
    val(seed)
   
  output:
    // tuple val(pool_id), path("${vireo_fixed_vcf}"), path("${vireo_fixed_vcf}.tbi"), emit: gt_pool
    path("${name}"),emit:output_collection
  script:
    file1 = files_to_concentrate[0]
    files= files_to_concentrate.join(" ")
    if (seed!='' ){
        seed_inset = "-seed '${seed}'"
    }else{
        seed_inset = ""
    }
    if (header_size>0){
        header_insert = "head -n ${header_size} ${file1} > ${name}"
        
    }else{
        header_insert = ""
    }

    header_insert

    """
        combine.py -d '${files}' -o ${name} ${seed_inset}
    """    
}