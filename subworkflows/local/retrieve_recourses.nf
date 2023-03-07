process RETRIEVE_RECOURSES{
  label 'process_tiny'
    // In nf there is a function collectFile - however if you provide a symlinked file directory tusing nf function will overwrite the source instead of replacing the file
    // This snipped is a replication of the function but as a nf module and hence the problem is avoided.
  publishDir "${params.outdir}",  mode: "${params.copy_mode}", overwrite: true  
  output:
    path("recourses"),emit:recourses
  script:
    if (params.reference_assembly_fasta_dir='"https://yascp.cog.sanger.ac.uk/public/10x_reference_assembly"'){
        get_genome = 'mkdir recourses/10x_reference_assembly && wget https://yascp.cog.sanger.ac.uk/public/10x_reference_assembly/genome.fa && wget https://yascp.cog.sanger.ac.uk/public/10x_reference_assembly/genome.fa.fai && mv genome.fa recourses/10x_reference_assembly/genome.fa && mv genome.fa.fai recourses/10x_reference_assembly/genome.fa.fai'
    }else{
        get_genome = ""
    }

    """
        mkdir recourses
        $get_genome
    """    
}