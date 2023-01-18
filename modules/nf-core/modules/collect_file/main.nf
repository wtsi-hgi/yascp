process collect_file{
  label 'process_tiny'
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
        seed_inset = "echo ${seed} > ${name}"
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
        ${seed_inset }
        ${header_insert}
        for OUTPUT in ${files}
        do
            echo ${files_to_concentrate}
            tail -n +${header_size+1} -q \$OUTPUT >> ${name}
        done
    """    
}