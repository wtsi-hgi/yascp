
process VIREO {
    tag "${samplename}"
    label 'process_high'
    publishDir "${params.outdir}/deconvolution/vireo/${samplename}/", mode: "${params.vireo.copy_mode}", overwrite: true,
	  saveAs: {filename -> filename.replaceFirst("vireo_${samplename}/","") }
    
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "/software/hgi/containers/mercury_scrna_deconvolution_62bd56a-2021-12-15-4d1ec9312485.sif"
        //// container "/software/hgi/containers/mercury_scrna_deconvolution_latest.img"
    } else {
        container "mercury/scrna_deconvolution:62bd56a"
    }

     when: 
      params.vireo.run

    input:
      tuple val(samplename), path(cell_data), val(n_pooled), path(donors_gt_vcf)
      

    output:
      tuple val(samplename), path("vireo_${samplename}/*"), emit: output_dir
      tuple val(samplename), path("vireo_${samplename}/donor_ids.tsv"), emit: sample_donor_ids 
      path("vireo_${samplename}/${samplename}.sample_summary.txt"), emit: sample_summary_tsv
      path("vireo_${samplename}/${samplename}__exp.sample_summary.txt"), emit: sample__exp_summary_tsv

    script:

      if (params.run_with_genotype_input){
        vcf = " -d ${donors_gt_vcf}"
      }else{
         vcf = ""
      }


    """

      umask 2 # make files group_writable

      vireo -c $cell_data -N $n_pooled -o vireo_${samplename} ${vcf} -t GT --randSeed 1 --nInit 200

      # add samplename to summary.tsv,
      # to then have Nextflow concat summary.tsv of all samples into a single file:

      cat vireo_${samplename}/summary.tsv | \\
        tail -n +2 | \\
        sed s\"/^/${samplename}\\t/\"g > vireo_${samplename}/${samplename}.sample_summary.txt

      cat vireo_${samplename}/summary.tsv | \\
        tail -n +2 | \\
        sed s\"/^/${samplename}__/\"g > vireo_${samplename}/${samplename}__exp.sample_summary.txt

    """
}
