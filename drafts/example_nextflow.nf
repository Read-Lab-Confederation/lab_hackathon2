nextflow.enable.dsl = 2
params.reads = '/home/ubuntu/nextflow/trial/input/*_R{1,2}.fastq.gz'
params.outdir = "./denovo"

reads = Channel.fromFilePairs(params.reads, checkIfExists: true)

// prints to the screen and to the log
log.info """
         Denovo Pipeline (version 1)
         ===================================
         reads        : ${params.reads}
         outdir       : ${params.outdir}
         """
         .stripIndent()


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

process assembly {
    /* 
       assembly step. here we used a conditional logic to choose the assembler
    */
    tag { sample_id }
    
    publishDir "$params.outdir/assemblies/", 
        mode: 'copy'
    
    input:
    tuple val(sample_id), path(filtered)  
    
    
    output:
    tuple val(sample_id), path("${sample_id}_contigs.fasta")

    script:
    """
    /home/ubuntu/spades/SPAdes-3.15.4-Linux/bin/spades.py -1 ${filtered[0]} \
    -2 ${filtered[1]} -o spades_output \
    -t 8 --isolate
    mv ./spades_output/contigs.fasta ${sample_id}_contigs.fasta

    """
}

workflow {
    fastp( reads )
    assembly( fastp.out.filtered )
}