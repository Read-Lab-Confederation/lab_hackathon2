nextflow.enable.dsl = 2

params.reads = '/home/ubuntu/extraVol/hackaton2/data/fastqs/*_R{1,2}.fastq.gz'
params.kraken_db = "/home/ubuntu/extraVol/hackaton2/data"
params.amr_fasta = "/home/ubuntu/extraVol/hackaton2/data/AMR_CDS"
params.data_dir = "/home/ubuntu/extraVol/hackaton2/data"
params.outdir = "/home/ubuntu/extraVol/hackaton2/output"
params.amr_dir = "/home/ubuntu/extraVol/hackaton2/output/AMR_CDS"


// prints to the screen and to the log
log.info """
         Denovo Pipeline (version 1)
         ===================================
         reads        : ${params.reads}
         outdir       : ${params.outdir}
         """
         .stripIndent()


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
      --detect_adapter_for_pe -w 8 -j ${sample_id}.fastp.json \
      --adapter_fasta /home/ubuntu/trimmomatic/Trimmomatic-0.39/adapters/NexteraPE-PE.fa

    """  
}  

process run_kma {
    tag "AMR on $sample_id"
    publishDir params.outdir, mode:'copy'

    input:
    tuple val(sample_id), path(filtered_reads)

    output:
    tuple val(sample_id), path("${sample_id}_kmamapped*")
    tuple val(sample_id), path("${sample_id}_kmamapped.res"), emit: kma_res
    tuple val(sample_id), path("${sample_id}_kmamapped.mapstat"), emit: kma_mapstat

    script:
    """
    
    kma -ipe  ${filtered_reads[0]}  ${filtered_reads[1]}  -t_db /home/ubuntu/extraVol/hackaton2/data/AMR_CDS -o ${sample_id}_kmamapped -mem_mode -ef -1t1 -cge -nf -vcf -t 1
    
    """
}

process format_kma_res {
  publishDir params.outdir, mode:'copy'
  
  input:
  tuple val(sample_id), path(kma_res)
  tuple val(sample_id), path(kma_mapstat)
  
  output:
  //path("${sample_id}.joined.tab")
  //tuple val(sample_id), path("${sample_id}.joined.name.tab")//, emit: kma_format
  path("${sample_id}.joined.name.tab")
  
  script:
  """
   awk '{if ((\$6 >= 20) && (\$7 >= 50)) print}'  ${kma_res} > ${sample_id}.res.filtered.tab
   tail -n +7 ${kma_mapstat} > ${sample_id}.mapstat.filtered.tab
   csvtk join -t -C '\$' -f 1 ${sample_id}.res.filtered.tab ${sample_id}.mapstat.filtered.tab > ${sample_id}_cur_joined.tab
   awk '{print FILENAME (NF?"\t":"") \$0}' *_cur_joined.tab | sed -E '2,\${/Template/d;}' > ${sample_id}.joined.tab
   awk -F'\t' -v value="${sample_id}" 'BEGIN{OFS=FS} {\$1 = value FS \$1; print}' ${sample_id}_cur_joined.tab > ${sample_id}.joined.name.tab
  
  """
}

process join_kma {
  
  publishDir params.outdir, mode:'copy'
  
  input:
  path '*.joined.tab'
  
  output:
  path 'all_kma_numerators_raw.tab'
  
  script:
  """
  awk '{print FILENAME (NF?"\t":"") \$0}' *.joined.tab | sed -E '2,\${/Template/d;}' > all_kma_numerators_raw.tab
  
  """
  
  
}

process run_kraken {
  
    publishDir params.outdir, mode:'copy'

    input:
    tuple val(sample_id), path(filtered_reads)
    
    output:
    path("${sample_id}_kraken_report_name.txt")

    
    script:
    """
    
    kraken2 --db  /home/ubuntu/extraVol/hackaton2/data/  --threads 8 --paired ${filtered_reads[0]} ${filtered_reads[1]} --output ${sample_id}_kraken_output.txt --report ${sample_id}_kraken_report.txt
    awk -F'\t' -v value="${sample_id}" 'BEGIN{OFS=FS} {\$1 = value FS \$1; print}' ${sample_id}_kraken_report.txt > ${sample_id}_kraken_report_name.txt
    
    
    """
  
}

process join_py {
  
  publishDir params.outdir, mode:'copy'
  
  input:
  tuple val(sample_id), path("*joined.tab")
  
  output:
  path "kma_comb.csv"
  
  
  script:
  """
  #!/home/ubuntu/miniconda3/bin/python
 
  import pandas as pd
 
  combined_table = pd.DataFrame()
  
  
  table_df = pd.read_csv("${kma_format}", sep = "\t")
  table_df.columns = table_df.columns.str.replace('#', '')
  table_df['Sample_id'] = "${sample_id}"
  combined_table = pd.concat([combined_table, table_df], ignore_index=True)
  combined_table.to_csv("kma_comb.csv", index = False)

  
  """
}

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
  echo The path is: $PATH
  Rscript --vanilla  ~/extraVol/hackaton2/programs/abundanceTable1.R
  
  """
  
}

workflow {
  reads = Channel.fromFilePairs(params.reads, checkIfExists: true)
  reference_fasta = Channel.fromPath(params.amr_fasta)
  index_amr( reference_fasta )
  filtered_ch = fastp( reads )
  run_kma(filtered_ch)
  // format_kma_res(run_kma.out.kma_res, run_kma.out.kma_mapstat) // | collect | join_kma
  format_kma_res(run_kma.out.kma_res, run_kma.out.kma_mapstat) | collect | join_cat
  //join_kma(format_kma_res.out.kma_join)
  run_kraken(filtered_ch) | collect | join_kraken
  abundance_tab(join_cat.out.kma_tab, join_kraken.out.kraken_tab)


  
}