#!/bin/bash

## This command calls the tool KMA to produced  the read counts associated with the gene hits
kma -ipe  ${fastq}_R1.fastq  ${fastq}_R2.fastq -t_db AMR_CDS_kmadb/AMR_CDS -o ${fastq}_kmamapped -mem_mode -ef -1t1 -cge -nf -vcf -t 1    


##Filters the results
 awk '{if (($6 >= 20) && ($7 >= 50)) print}'  ${sample}_KMA_mapped.res > ${sample}_KMA_mapped.res.filtered.tab

## filters the mapstat file
tail -n +7 ${sample}.mapstat > ${sample}.mapstat.filtered.tab

## copies the mapstat and res files

awk '{print FILENAME (NF?"\t":"") $0}' *kmamapped.joined.tab | sed -E '2,${/Template/d;}' > all_kma_numerators_raw.tab

## command to run krakenn (assuming you already have a database
kraken2 --db  kraken_db/  --threads 4 --paired ${fastq}_R1.fastq ${fastq} _R2.fastq --output kraken2_results/${fastq}_output.txt --report kraken2_results/${fastq}_kraken_report

## prints the file name to the first column and then concatenates all of the report files together
Awk '{print FILENAME (NF?"\t":"") $0}' *kraken_report > all_kraken_reports.tab

## filters only those rows where the taxonomic classification was bacteria (2) or unclassified (0)
awk '{if (($6 == 2) || ($6 == 0)) {print} }' all_kraken_reports.tab > all_kraken_bacteria.tab

## congrats now you can work with the R code to actually calculate some stuff
