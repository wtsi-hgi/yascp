process SUBSET_BAM_PER_BARCODES_AND_VARIANTS {
    tag "${samplename}"    
    label 'process_low'
   
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "/lustre/scratch123/hgi/projects/ukbb_scrna/pipelines/Pilot_UKB/carls_data/analysis/bam_tool_processing_05_04_2024.sif"
    } else {
        container " mercury/bam_tool_processing:05_04_2024"
    }
    
    input:
        tuple val(sample), path(barcodes), path(vcf), path(bam)

    output:
        tuple val(sample), val("${donor}") ,path("${sample}_filtered.bam"),path("${sample}_filtered.bam.csi"), path(barcodes), emit: freebayes_input

    script:
        def donor_split = "${barcodes}".tokenize('.')
        donor = donor_split[0]
        """ 
            echo ${donor}
            bcftools sort ${vcf} -Oz -o srt_${vcf}
            filter_bam_file_for_popscle_dsc_pileup.sh ${bam} ${barcodes} srt_${vcf} ${sample}_filtered.bam
        """

}