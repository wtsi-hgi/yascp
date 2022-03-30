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
      tuple val(pool_id), path("${cellranger_possorted_bam}"), path("${outdir}"), emit: poolid_bam_dirpath
      path("${oufofn}", emit: bcfiles_fofn)

    script:
    oufnprfx="${pool_id}_barcodes"
    outdir="${pool_id}_bcfiles"
    oufofn="${pool_id}_fofn.txt"
    """
      write_barcode_to_file () {
          barcode=\$1
          donor=\$2
          echo "\${barcode}" >> ${outdir}/${oufnprfx}_\${donor}.txt
      }

      rm -fr ${outdir}
      mkdir ${outdir}
      cut -f1,2 ${vireo_donor_barcode_tsv} | tail -n +2 |\
      while read line
      do
        write_barcode_to_file \${line}
      done

      echo "pool_id,barcode_file" > ${oufofn}
      for fn in \$(ls ${outdir}/${oufnprfx}_*.txt)
      do
        fnam=\$(basename \${fn})
        echo "${pool_id},\${fnam}" >> ${oufofn}
      done
    """
}

process ADD_PATH_TO_FOFN
{
  label 'process_tiny'

  input:
    tuple val(pool_id),

  output:
    path ("${updated_fofn}", emit:full_path_fofn)

  script:
  updated_fofn = "${bcfiles_fofn}".minus(".txt").plus("_fullpaths.txt")
  """
  dirpath=\$(dirname ${bcfiles_fofn})
  awk --assign dp=\${dirpath} -F"," '{if(NR>1){printf"%s,%s/%s,%s/%s\\n",\$1,dp,\$2,dp,\$3}else{print \$0}}' ${bcfiles_fofn} > ${updated_fofn}
  """
}

process SPLIT_BAM_BY_CELL_BARCODES
{
    //label 'process_medium'

    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
      container "/software/hgi/containers/wtsihgi-nf_yascp_htstools-1.0.sif"
    } else {
      container "mercury/wtsihgi-nf_yascp_htstools-1.0"
    }

    publishDir path: "${params.outdir}/bam/${pool_id}",
               mode: "${params.copy_mode}",
               overwrite: "true"

    input:
      tuple val(pool_id), path(cellranger_possorted_bam), path(dirpath), val(vireo_donor_barcode_filnam)

    output:
      path("${oufnprfx}_possorted_bam.bam", emit: possorted_bam_files)

    script:
    oufnprfx = "${vireo_donor_barcode_filnam}".minus(".txt")
    bcfilpath = "${dirpath}/${vireo_donor_barcode_filnam}"
    """
      echo "pool_id = ${pool_id}"
      echo "cellranger_possorted_bam = ${cellranger_possorted_bam}"
      echo "dirpath = ${dirpath}"
      echo "vireo_donor_barcode_filnam = ${vireo_donor_barcode_filnam}"
      samtools view --threads ${task.cpus} --tag-file CB:${bcfilpath} --bam -o ${oufnprfx}_possorted_bam.bam ${cellranger_possorted_bam}
    """
}

workflow split_bam_by_donor
{
  take:
    poolid_bam_vireo_donor_barcodes

  main:
    SPLIT_CELL_BARCODES_PER_DONOR(poolid_bam_vireo_donor_barcodes)

    SPLIT_CELL_BARCODES_PER_DONOR.out.bcfiles_fofn
    .splitCsv(header: true)
    .map { row -> [row.pool_id, row.barcode_file]}
    .set {ch_bcfiles_fofn}

    SPLIT_CELL_BARCODES_PER_DONOR.out.poolid_bam_dirpath
    .set {ch_poolid_bam_dirpath}

    ch_bcfiles_fofn.subscribe onNext: { println "ch_bcfiles_fofn = ${it}" }, onComplete: { println "Done.\n" }
    ch_poolid_bam_dirpath.subscribe onNext: {"println ch_poolid_bam_dirpath = ${it}" }, onComplete: { println "Done.\n" }

    ch_poolid_bam_dirpath.combine(ch_bcfiles_fofn, by: 0)
    .set {ch_poolid_bam_dirpath_bcfile}

    ch_poolid_bam_dirpath_bcfile
    .subscribe onNext: { println "ch_poolid_bam_dirpath_bcfile = ${it}" }, onComplete: { println "Done.\n" }

    SPLIT_BAM_BY_CELL_BARCODES(
      ch_poolid_bam_dirpath_bcfile
    )

  emit:
    possorted_bam_files = SPLIT_BAM_BY_CELL_BARCODES.out.possorted_bam_files

}
