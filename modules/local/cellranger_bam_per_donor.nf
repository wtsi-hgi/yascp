process SPLIT_CELL_BARCODES_PER_DONOR
{
    label 'process_tiny'

    publishDir path: "${params.outdir}/bam/${pool_id}",
               mode: "${params.copy_mode}",
               pattern: "${oufnprfx}_*.txt",
               overwrite: "true"

    input:
      tuple val(pool_id), path(cellranger_possorted_bam), path(vireo_donor_barcode_tsv)

    output:
      path "${oufnprfx}_*.txt", emit: per_donor_barcode_files
      tuple val(pool_id), path(cellranger_possorted_bam), emit: experiment_possorted_bam

    script:
    oufnprfx="${pool_id}_barcodes"
    """
      write_barcode_to_file () {
          barcode=\$1
          donor=\$2
          echo "\${barcode}" >> ${oufnprfx}_\${donor}.txt
      }

      rm -f ${oufnprfx}_*.txt

      cut -f1,2 ${vireo_donor_barcode_tsv} | tail -n +2 |\
        while read line
      do
        write_barcode_to_file \${line}
      done

    """
}

process SPLIT_BAM_BY_CELL_BARCODES
{
    label 'process_low'

    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
      container "/software/hgi/containers/wtsihgi-nf_yascp_htstools-1.0.sif"
    } else {
      container "mercury/wtsihgi-nf_yascp_htstools-1.0"
    }

    publishDir path: "${params.outdir}/bam/${pool_id}",
               mode: "${params.copy_mode}",
               overwrite: "true"

    input:
      tuple val(pool_id), path(cellranger_possorted_bam)
      path(vireo_donor_barcode_file)

    output:
      path "${oufnprfx}_possorted_bam.bam", emit: possorted_bam_files

    script:
    oufnprfx = "${vireo_donor_barcode_file}".minus(".txt")
    """
      samtools view --tag-file CB:${vireo_donor_barcode_file} --bam -o ${oufnprfx}_possorted_bam.bam ${cellranger_possorted_bam}
    """
}

workflow SPLIT_BAM_PER_DONOR
{
  take:
    poolid_bam_vireo_donor_barcodes

  main:
    SPLIT_CELL_BARCODES_PER_DONOR(poolid_bam_vireo_donor_barcodes)
    SPLIT_CELL_BARCODES_PER_DONOR.out.per_donor_barcode_files
      .flatten()
      .set { ch_per_donor_barcode_files }

    SPLIT_BAM_BY_CELL_BARCODES(
      SPLIT_CELL_BARCODES_PER_DONOR.out.experiment_possorted_bam,
      ch_per_donor_barcode_files
      )

  emit:
    possorted_bam_files = SPLIT_BAM_BY_CELL_BARCODES.out.possorted_bam_files
}
