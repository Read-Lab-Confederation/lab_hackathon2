### Set up the R script ###
## Input tables be recognized by R 

library(optparse)
library(tidyr)
library(dplyr)
library(stringr)
library(ggplot2)

option_list <- list( 
  make_option(c("-i", "--input"),  type="character", metavar = "FILE", ), 
  make_option(c("-d", "--input_denominator"),  type="character", metavar = "FILE", ))
opt=parse_args(OptionParser(option_list=option_list))

numerator<-opt$input
demoninator<-opt$input_denominator

KMA_tables<- read.table(numerator, sep = "\t", fill = TRUE)
bacteria_kraken_table<- read.table(demoninator, sep = "\t", fill = TRUE)


## Functions to make R script recognized by bash 


## For KMA specifically, samples were removed if the filtered "res" file or the "mapstat" file for the same sequence was blank
##

## KMA Can give you the 
KMA_tables <- read.table("~/tiramisu/github_work/lab_hackathon2/data/all_kma_numerators_raw.tab", header = TRUE, comment.char = "$", sep ='\t', fill = TRUE)

colnames(KMA_tables)[1] <- "SAMPLE_ID"

## FILTER OUT THE HEADERS, the excluded sequences  ##
KMA_tables$SAMPLE_ID <- str_remove(KMA_tables$SAMPLE_ID,"_.*")
KMA_tables <- KMA_tables %>% filter(!(X.Template %in% c("#Template","# refSequence", "Template")))

## TO CALCULATE THE AMR RPKM ##

KMA_tables <- KMA_tables %>%  mutate(template_length_kb = as.numeric(Template_length)/1000)

## NOW NEED TO READ IN THE BACTERIAL READS COUNT ## 

bacteria_kraken_table <- read.table("~/tiramisu/github_work/lab_hackathon2/data/kraken2_results/all_kraken_bacteria.tab", header = FALSE, comment.char = "$", sep ='\t', fill = TRUE)

colnames(bacteria_kraken_table)[1] <- "SAMPLE_ID"
colnames(bacteria_kraken_table)[2] <- "percent_rooted_reads"
colnames(bacteria_kraken_table)[3] <- "number_rooted_reads"
colnames(bacteria_kraken_table)[4] <- "number_taxon_reads"
colnames(bacteria_kraken_table)[5] <- "taxon_rank"
colnames(bacteria_kraken_table)[6] <- "taxon_symbol"
colnames(bacteria_kraken_table)[7] <- "taxon"

bacteria_kraken_table$SAMPLE_ID <- str_remove(bacteria_kraken_table$SAMPLE_ID,"_.*")
bacteria_kraken_table2 <- bacteria_kraken_table %>% filter(taxon_symbol == 2)


abundance_table <- left_join(KMA_tables, bacteria_kraken_table2, by = "SAMPLE_ID")

abundance_table$RPKM <- abundance_table$readCount/(abundance_table$template_length_kb*abundance_table$number_rooted_reads)*1000000000

abundance_table <- abundance_table %>% separate(X.Template, into = c("x1", "x2", "x3", "num1","num2","Gene_Symbol","Gene_symbol2","description"), sep = "\\|")

write.table(abundance_table, file=paste("gene_abundance_table.tab", sep=""), sep = "\t", quote = FALSE)

abundance_plot <- ggplot(abundance_table, aes(x= reorder(Gene_Symbol, -RPKM), y= RPKM)) +
  geom_bar(stat = "identity")

ggsave("abundance_plot.png", plot = abundance_plot)