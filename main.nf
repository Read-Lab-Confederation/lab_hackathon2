#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

// set paths
params.reads = "$baseDir/data/*_R{1,2}.fastq.gz"
params.amr_fasta = "$baseDir/data/AMR_CDS"
params.data_dir = "$baseDir/data"
params.outdir = "$baseDir/output"
params.help = ""

def helpMessage() {
  log.info """
        Add Help Menu!!
        """ 
}

// Show help message
if (params.help) {
    helpMessage()
    exit 0
}


// prints to the screen and to the log
log.info """
         Denovo Pipeline (version 1)
         ===================================
         reads        : ${params.reads}
         outdir       : ${params.outdir}
         """
         .stripIndent()

// index amr fasta
process index_amr {

    publishDir params.data_dir, mode:'copy'
    
    input:
    path fasta_file
    
    output:
    file "*"
    
    script:
    
    """
    kma index -i ${fasta_file} -o AMR_CDS
    
    """
  
}

// filter reads with basic fastp
process fastp {
    /* 
       fastp process to remove adapters and low quality sequences
    */
    tag "filter $sample_id"

    input:
    tuple val(sample_id), path(reads) 
    
    output:
    tuple val(sample_id), path("${sample_id}_filt_R{1,2}.fastq.gz")

 
    script:
    """
    fastp -i ${reads[0]} -I ${reads[1]} \
      -o ${sample_id}_filt_R1.fastq.gz -O ${sample_id}_filt_R2.fastq.gz \
      --detect_adapter_for_pe -w 8 -j ${sample_id}.fastp.json

    """  
}  

// run kma 
process run_kma {
    tag "AMR on $sample_id"

    input:
    tuple val(sample_id), path(filtered_reads)

    output:
    tuple val(sample_id), path("${sample_id}_kmamapped.res"), emit: kma_res
    tuple val(sample_id), path("${sample_id}_kmamapped.mapstat"), emit: kma_mapstat

    script:
    """
    
    kma -ipe  ${filtered_reads[0]}  ${filtered_reads[1]}  -t_db $baseDir/data/AMR_CDS -o ${sample_id}_kmamapped -mem_mode -ef -1t1 -cge -nf -vcf -t 1
    
    """
}

// join kma results with csvtk and add sample name
process format_kma_res {
  
  input:
  tuple val(sample_id), path(kma_res)
  tuple val(sample_id), path(kma_mapstat)
  
  output:
  
  path("${sample_id}.joined.name.tab")
  
  script:
  """
   awk '{if ((\$6 >= 20) && (\$7 >= 50)) print}'  ${kma_res} > ${sample_id}.res.filtered.tab
   tail -n +7 ${kma_mapstat} > ${sample_id}.mapstat.filtered.tab
   csvtk join -t -C '\$' -f 1 ${sample_id}.res.filtered.tab ${sample_id}.mapstat.filtered.tab > ${sample_id}_cur_joined.tab
   awk -F'\t' -v value="${sample_id}" 'BEGIN{OFS=FS} {\$1 = value FS \$1; print}' ${sample_id}_cur_joined.tab > ${sample_id}.joined.name.tab
  
  """
}

// run kraken2
process run_kraken {

    input:
    tuple val(sample_id), path(filtered_reads)
    
    output:
    path("${sample_id}_kraken_report_name.txt")

    
    script:
    """
    
    kraken2 --db  $baseDir/data/ --threads 8 --paired ${filtered_reads[0]} ${filtered_reads[1]} --output ${sample_id}_kraken_output.txt --report ${sample_id}_kraken_report.txt
    awk -F'\t' -v value="${sample_id}" 'BEGIN{OFS=FS} {\$1 = value FS \$1; print}' ${sample_id}_kraken_report.txt > ${sample_id}_kraken_report_name.txt
    
    
    """
  
}

// join formatted kma tables for all samples
process join_cat {
  publishDir params.outdir, mode:'copy'
  input:
  path "*.joined.name.tab"
  output:
  path 'kma_all_joinned.tsv', emit: kma_tab
  script:
  """
   cat *.joined.name.tab > kma_all_joinned.tsv
  
  """
  
}

// combine all samples of kraken result tables 
process join_kraken {
  publishDir params.outdir, mode:'copy'
  input:
  path "*_name.txt"
  output:
  path 'kraken_all_report.tsv', emit: kraken_tab
  script:
  """
   cat *_name.txt > kraken_all_report.tsv
  
  """
  
}

// format abundance table using R program
process abundance_tab {
  publishDir params.outdir, mode:'copy'
  
  input:
  path kma_tab
  path kraken_tab
  
  output:
  path 'gene_abundance_table.tsv'
  path 'abundance_plot.png'
  
  script:
  """

  abundanceTable1.R
  
  """
  
}

// workflow 
workflow {
  reads = Channel.fromFilePairs(params.reads, checkIfExists: true)
  reference_fasta = Channel.fromPath(params.amr_fasta)
  index_amr( reference_fasta )
  filtered_ch = fastp( reads )
  run_kma(filtered_ch)
  format_kma_res(run_kma.out.kma_res, run_kma.out.kma_mapstat) | collect | join_cat
  run_kraken(filtered_ch) | collect | join_kraken
  abundance_tab(join_cat.out.kma_tab, join_kraken.out.kraken_tab)

}
