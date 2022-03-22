process SPLIT_BAM_PER_DONOR
{
    label 'process_tiny'

    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
      container "/software/hgi/containers/wtsihgi-nf_yascp_htstools-1.0.sif"
    } else {
      container "mercury/wtsihgi-nf_yascp_htstools-1.0"
    }

    input:
      path(vireo_donor_barcode_tsv)
      path(cellranger_possorted_bam)

    script:
    """
      write_barcode_to_file () {
          barcode=\$1
          donor=\$2
          echo "\${barcode}" >> barcodes_\${donor}.txt
      }

      rm -f barcodes_*.txt

      cut -f1,2 ${vireo_donor_barcode_tsv} | tail -n +2 |\
        while read line
      do
        write_barcode_to_file \${line}
      done

      for bcfn in \$(ls barcodes_*.txt)
      do
          oufn=\$(basename \${bcfn} .txt)
          echo "output file: \${oufn}"
          samtools view --tag-file CB:\${oufn}.txt --bam -o \${oufn}_possorted_bam.bam ${cellranger_possorted_bam}
      done

    """
}
