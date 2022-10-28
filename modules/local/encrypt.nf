process ENCRYPT_DIR
{
  label 'process_tiny'

  if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
    container "/software/hgi/containers/wtsihgi-nf_ega_cryptor-1.0.sif"
  } else {
    container "wtsihgi/egacryptor:fec9c55377f6" // old container -> "mercury/wtsihgi-nf_ega_cryptor-1.0"
  }

  publishDir path: "${params.outdir}/ukbb/",
             saveAs: {filnam -> file(filnam).getName()},
             mode: 'copy',
             overwrite: 'true'
  publishDir "${versionsDir}", pattern: "*.versions.yml", mode: "${params.versions.copy_mode}"

  input:
    path(input_dir)

  output:
    path("${outdir}/*", emit:encrypted_files)
    path("${outdir}", emit:encrypted_dir)
    // path("${outdir}/*/*.md5", emit:checksums)
    path ('*.versions.yml')         , emit: versions 

  script:
  outdir="encryptor_output"
  """
    version = bash encript.sh ${input_dir} ${outdir}
    
    ####
    ## capture software version
    ####
    echo "${task.process}:" > ${task.process}.versions.yml
    echo "    EGA-Cryptor: \$version" >> ${task.process}.versions.yml
  """
}
