#!/usr/bin/env Rscript 

library(tidyr)
library(dplyr)
library(stringr)
library(ggplot2)

print(list.files())

## KMA Can give you the 
KMA_tables <- read.table("kma_all_joinned.tsv", header = TRUE, comment.char = "$", sep ='\t', fill = TRUE)

colnames(KMA_tables)[1] <- "SAMPLE_ID"

## FILTER OUT THE HEADERS, the excluded sequences  ##
#KMA_tables$SAMPLE_ID <- str_remove(KMA_tables$SAMPLE_ID,"_.*")

KMA_tables <- KMA_tables %>% filter(!(X.Template %in% c("#Template","# refSequence", "Template")))

## TO CALCULATE THE AMR RPKM ##
KMA_tables <- KMA_tables %>%  mutate(template_length_kb = as.numeric(Template_length)/1000)

## NOW NEED TO READ IN THE BACTERIAL READS COUNT ## 

bacteria_kraken_table <- read.table("kraken_all_report.tsv", header = FALSE, comment.char = "$", sep ='\t', fill = TRUE)

colnames(bacteria_kraken_table)[1] <- "SAMPLE_ID"
colnames(bacteria_kraken_table)[2] <- "percent_rooted_reads"
colnames(bacteria_kraken_table)[3] <- "number_rooted_reads"
colnames(bacteria_kraken_table)[4] <- "number_taxon_reads"
colnames(bacteria_kraken_table)[5] <- "taxon_rank"
colnames(bacteria_kraken_table)[6] <- "taxon_symbol"
colnames(bacteria_kraken_table)[7] <- "taxon"

bacteria_kraken_table2 <- bacteria_kraken_table %>% filter(taxon_symbol == 2)


abundance_table <- left_join(KMA_tables, bacteria_kraken_table2, by = "SAMPLE_ID")

# make columns numric
abundance_table$readCount = as.numeric(as.character(abundance_table$readCount))
abundance_table$template_length_kb = as.numeric(as.character(abundance_table$template_length_kb))
abundance_table$number_rooted_reads = as.numeric(as.character(abundance_table$number_rooted_reads))

abundance_table$RPKM <- abundance_table$readCount/(abundance_table$template_length_kb*abundance_table$number_rooted_reads)*1000000000

abundance_table <- abundance_table %>% separate(X.Template, into = c("x1", "x2", "x3", "num1","num2","Gene_Symbol","Gene_symbol2","description"), sep = "\\|")

write.table(abundance_table, file="gene_abundance_table.tsv", sep = "\t", quote = T, row.names = F)

abundance_plot <- ggplot(abundance_table, aes(x= reorder(Gene_Symbol, -RPKM), y= RPKM)) +
  geom_bar(stat = "identity")

ggsave("abundance_plot.png", plot = abundance_plot)