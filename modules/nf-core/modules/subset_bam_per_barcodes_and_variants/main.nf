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


process SUBSET_BAM_PER_BARCODES{


    tag "${samplename}"    
    label 'process_low'
   
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "/lustre/scratch123/hgi/projects/ukbb_scrna/pipelines/Pilot_UKB/carls_data/analysis/bam_tool_processing_05_04_2024.sif"
    } else {
        container "mercury/bam_tool_processing:05_04_2024"
    }

    publishDir  path: "${params.outdir}/handover/Donor_Quantification/${sample}",
                mode: "${params.copy_mode}",
                overwrite: "true"

    input:

    // CRD_CMB12979963, /lustre/scratch123/hgi/mdt2/teams/hgi/mo11/tmp_projects/cardinal/qc_F3/UKBB_ELGH_5th_July_2022/work/5c/dada93d9cbebef97d30e4014e131cf/input__CRD_CMB12979963/possorted_genome_bam.bam, /lustre/scratch123/hgi/mdt2/teams/hgi/mo11/tmp_projects/cardinal/qc_F3/UKBB_ELGH_5th_July_2022/work/5c/dada93d9cbebef97d30e4014e131cf/input__CRD_CMB12979963/possorted_genome_bam.bam.bai, /lustre/scratch123/hgi/mdt2/teams/hgi/mo11/tmp_projects/cardinal/qc_F3/UKBB_ELGH_5th_July_2022/results/nf-preprocessing/cellbender/CRD_CMB12979963/cellbender-epochs_250__learnrt_0pt000005__zdim_100__zlayer_500__lowcount_10/cellbender-FPR_0pt1-filtered_10x_mtx/barcodes.tsv.gz
        tuple val(sample), val(sample_donor),path(barcodes), path(bam), path(bai), path(cellbender_barcodes)
        path(genome)
    output:
        tuple val(sample),path("${sample_donor}.cram"), emit: cram_files
    //     tuple val(sample), val("${donor}") ,path("${sample}_filtered.bam"),path("${sample}_filtered.bam.csi"), path(barcodes), emit: freebayes_input

    script:

        """ 
            samtools view --threads ${task.cpus} --tag-file CB:${barcodes} \
                --cram -T ${genome}/genome.fa \
                -o ${sample_donor}.cram ${bam}
        """

}