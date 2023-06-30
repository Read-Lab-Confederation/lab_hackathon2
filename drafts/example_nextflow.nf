nextflow.enable.dsl = 2

params.reads = '/home/ubuntu/extraVol/hackaton2/data/fastqs/*_R{1,2}.fastq.gz'
params.kraken_db = "/home/ubuntu/extraVol/hackaton2/data"
params.amr_fasta = "/home/ubuntu/extraVol/hackaton2/data/AMR_CDS"
params.outdir = "/home/ubuntu/extraVol/hackaton2/data"



// prints to the screen and to the log
log.info """
         Denovo Pipeline (version 1)
         ===================================
         reads        : ${params.reads}
         outdir       : ${params.outdir}
         """
         .stripIndent()


process index_amr {
      /* 
       index AMR fasta file
    */
    
    publishDir "$params.outdir/", 
        mode: 'copy'
    /*
    conda: 'kma'
    */
    
    input:
    path fasta_file
    
    output:
    path AMR_CDS, emit: index_files
    
    script:
    
    """
    kma index -i ${fasta_file} -o AMR_CDS
    
    """
  
}


process fastp {
    /* 
       fastp process to remove adapters and low quality sequences
    */
    tag "filter $sample_id"

    input:
    tuple val(sample_id), path(reads) 
    
    output:
    tuple val(sample_id), path("${sample_id}_filt_R{1,2}.fastq.gz"), emit: filtered
    path("${sample_id}.fastp.json"), emit: json

 
    script:
    """
    fastp -i ${reads[0]} -I ${reads[1]} \
      -o ${sample_id}_filt_R1.fastq.gz -O ${sample_id}_filt_R2.fastq.gz \
      --detect_adapter_for_pe -w 8 -j ${sample_id}.fastp.json \
      --adapter_fasta /home/ubuntu/trimmomatic/Trimmomatic-0.39/adapters/NexteraPE-PE.fa

    """  
}  


process run_kma {
  
  
  input:
  tuple val(sample_id), path(filtered), path(index_files)
  
  
}


workflow {
  input_channel = Channel.fromPath(params.amr_fasta)
  reads = Channel.fromFilePairs(params.reads, checkIfExists: true)
  index_amr( input_channel )
 
  fastp( reads )
  /*  
    assembly( fastp.out.filtered )
  */
}