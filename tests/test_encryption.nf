#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

// include { ENCRYPT_DIR } from '../modules/local/encrypt.nf'

process ENCRYPT_DIR
{
  label 'process_tiny'

  if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
    container "/software/hgi/containers/wtsihgi-nf_ega_cryptor-1.0.sif"
  } else {
    container "mercury/wtsihgi-nf_ega_cryptor-1.0"
  }

  publishDir path: "${params.outdir}/ukbb/",
             saveAs: { filnam -> file(filnam).getName() },
             mode: 'copy',
             overwrite: "true"

  input:
    path(input_dir)

  output:
    path("${outdir}/*.gpg", emit:outfil)
    path("${outdir}/*.md5", emit:checksum)

  script:
  outdir="encryptor_output"
  """
    mkdir ${outdir}
    echo "test1.gpg" > ${outdir}/test1.gpg
    echo "test1.md5" > ${outdir}/test1.md5
  """
}

workflow TEST_ENCRYPT_DIR {
  Channel.fromPath("/lustre/scratch123/hgi/projects/ukbb_scrna/pipelines/hp3_dev/elgh_yascp/test/")
  .set {ch_dir}
  ENCRYPT_DIR(ch_dir)
}
