#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

// set paths
params.reads = "$baseDir/data/*_R{1,2}.fastq.gz"
params.amr_fasta = "$baseDir/data/compiled_AMR_database.fasta"
params.data_dir = "$baseDir/data/index"
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

    //publishDir params.data_dir, mode:'copy'
    
    input:
    path fasta_file
    
    output:
    file "*"
    
    script:
    //-o compiled_AMR
    """
    kma index -i ${fasta_file} -o compiled_AMR
    
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
    
    kma -ipe  ${filtered_reads[0]}  ${filtered_reads[1]}  -t_db $baseDir/data/index/compiled_AMR -o ${sample_id}_kmamapped -mem_mode -ef -1t1 -cge -nf -vcf -t 1
    
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
   awk '{if ((\$6 >= 90) && (\$7 >= 60)) print}'  ${kma_res} > ${sample_id}.res.filtered.tab
   tail -n +7 ${kma_mapstat} > ${sample_id}.mapstat.filtered.tab
   csvtk join -t -C '\$' -f 1 ${sample_id}.res.filtered.tab ${sample_id}.mapstat.filtered.tab > ${sample_id}_cur_joined.tab
   awk -F'\t' -v value="${sample_id}" 'BEGIN{OFS=FS} {\$1 = value FS \$1; print}' ${sample_id}_cur_joined.tab > ${sample_id}.joined.name.tab
  
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

// run kraken2
process run_kraken {
    //publishDir params.outdir, mode:'copy'
    
    input:
    tuple val(sample_id), path(filtered_reads)
    
    output:
    tuple val(sample_id), path("${sample_id}.kr"), emit: kraken_report
    
    script:
    """
    kraken2 --db  $baseDir/data/ --threads 8 --paired ${filtered_reads[0]} ${filtered_reads[1]} --output ${sample_id}_kraken_output.txt --report ${sample_id}.kr
    """  
}

// edit kraken out
process edit_kraken{
    //publishDir params.outdir, mode:'copy'

    input:
    tuple val(sample_id), path(kraken_report)

    output:
    path("${sample_id}_kraken_report_name.txt")

    script:
    """
    awk -F'\t' -v value="${sample_id}" 'BEGIN{OFS=FS} {\$1 = value FS \$1; print}' ${kraken_report} > ${sample_id}_kraken_report_name.txt
    """
}

// run Bracken
process run_bracken {
  publishDir params.outdir, mode:'copy'

  input:
  tuple val(sample_id), path(kraken_report)

  output:
  path("${sample_id}.br")

  script:
  """
  $baseDir/Bracken/bracken -d $baseDir/data/ -i ${kraken_report} -o ${sample_id}.br -r 150 -l "S"
  """
}

// create manifest for bracken files
process create_bracken_manifest {
  //publishDir params.outdir, mode:'copy'
 
  input:
  //path "*.br"
  path params.outdir

  output:
  path 'bracken_manifest.tsv', emit: bracken_mani

  script:
  """
    $baseDir/scripts/bracken_mnfst_gen.sh ${params.outdir}
  """
} 

// create pdf of the Fractional Read Abundance
process create_Kraken_Bracken_plot {
  publishDir params.outdir, mode:'copy'
 
  input:
  path(bracken_mani)

  output:
  path 'Fractional_Read_Abundance.pdf'
  path 'FRA_table.csv', emit: fra_table

  script:
  """
    python $baseDir/scripts/Kraken-Bracken-plot.py -i ${bracken_mani} -p Fractional_Read_Abundance -o FRA_table.csv
  """
} 

// create pdf of the FRA without the dominant species
process create_wo_dominant_plot {
  publishDir params.outdir, mode:'copy'

  input:
  path(fra_table)

  output:
  path 'partial_fractional_abundance.pdf'

  script:
  """
    python $baseDir/scripts/deeper_FRA_vis.py ${fra_table}
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
  path 'gene_abundance_table.tsv', emit: gene_abundance_table
  //path 'abundance_plot.png'
  
  script:
  """

  abundanceTable1.R
  
  """
  
}

// create vis of gene RPKM
process create_gene_RPKM_plot {
  publishDir params.outdir, mode:'copy'

  input:
  path(gene_abundance_table)

  output:
  path 'genes_rpkm.pdf'

  script:
  """
    python $baseDir/scripts/rpkm_vis.py ${gene_abundance_table}
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
  kraken_out_ch = run_kraken(filtered_ch)
  edit_kraken(kraken_out_ch) | collect | join_kraken
  abundance_tab(join_cat.out.kma_tab, join_kraken.out.kraken_tab)
  create_gene_RPKM_plot(abundance_tab.out.gene_abundance_table)
  filtered_brak_ch = run_bracken(kraken_out_ch) | collect | create_bracken_manifest
  create_Kraken_Bracken_plot(filtered_brak_ch)
  create_wo_dominant_plot(create_Kraken_Bracken_plot.out.fra_table)
}