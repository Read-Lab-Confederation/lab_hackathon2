# lab_hackathon2
Our project is to adapt Brooke's script for determining AMR reads per million from meta genomic data

Here is the [the google doc](https://docs.google.com/document/d/1a1NjFz8dDE0VPHwXtsbpe8BlrzB2dXce/edit) with project info.

[Lab meeting presentation](https://docs.google.com/presentation/d/1DkSnNEFyrNsgcvd66kASLn81-8fmhp7KXd5w81L1w_0/edit#slide=id.g17b97a35150_1_0) where we did git tutorial.

[Github site](https://github.com/Read-Lab-Confederation/github-collab-practice) for practising github pushes and pulls.

## January 2024 Hackathon

*Preparing the local environment*

** Note: as of 2024-01-30 there is a bug in the nextflow script when only ONE pair of fastqs is in the input directory**≈y

```git clone git@github.com:Read-Lab-Confederation/lab_hackathon2.git```

```conda create -c bioconda -n hack2 nextflow kma kraken2 csvtk fastp```

```conda activate hack2```

*add data (create data directory if it doesnt already exist)*

```cd lab_hackathon2/data/```

```wget https://ftp.ncbi.nlm.nih.gov/pathogen/Antimicrobial_resistance/AMRFinderPlus/database/latest/AMR_CDS```

```mkdir fastqs```

```cd fastqs```

```wget -O simulated_metagenome_R1.fastq.gz https://zenodo.org/record/6543357/files/simulated_metagenome_1.fq.gz?download=1```

```wget -O simulated_metagenome_R2.fastq.gz https://zenodo.org/record/6543357/files/simulated_metagenome_2.fq.gz?download=1```

```# cp small fastq files from Michael David metagenome project

cp /mnt/tiramisu/emergent/projects/SEMAPHORE/data/fastqs/semaphore/microbiome/data_files/S.190905.00152_* ./```

```cd ../```

```wget https://genome-idx.s3.amazonaws.com/kraken/k2_standard_08gb_20230605.tar.gz```

```tar -xvzf k2_standard_08gb_20230605.tar.gz```

## Agenda for thurs 29th June 2023

1. Clone the github repo
2. Create conda environments based on YAML. (Any other softwhere we need to add to the environment , like nextflow?)
   
```conda create -c bioconda -n hack2 nextflow kma kraken2 csvtk```
   
```conda activate hack2```
   
4. Download test data sets

```wget https://zenodo.org/record/6543357/files/simulated_metagenome_1.fq.gz?download=1```

```wget https://zenodo.org/record/6543357/files/simulated_metagenome_2.fq.gz?download=1```

5. Download kraken database and the AMR gene database

```wget https://genome-idx.s3.amazonaws.com/kraken/k2_standard_08gb_20230605.tar.gz```

(move the kraken database outside of your github directory)

```wget https://ftp.ncbi.nlm.nih.gov/pathogen/Antimicrobial_resistance/AMRFinderPlus/database/latest/AMR_CDS```

Create kma index

```kma index -i AMR_CDS -o templates```

##Agenda for friday 29th June 2023

1. Downsample the synthetic datasets?
2. One group work on Rscript
```Rscript abundanceTable.r --input all_kma_numerators_raw.tab --input_denominator all_kraken_bacteria.tab```
3. The other group work on converting the existing shell script to nextflow
4. Docker container?


##Andrei's work on the pipeline
  1. Results and instructions:
  
    a) wrapped the whole pipeline in nextflow
    
    b) created new yaml for conda dependencies
    
    c) to install the environment type: conda env create -fhack2_nextflow.yaml
    
    d) to run the pipeline type:  bash -i runAll.sh
    
    e) current version uses local paths for references, so data/, bin/, main.nf, and runAll.sh should be in the same directory
    
    f) data/ directory can be downloaded from s3://transfer-files-emory/amrKma/data.tar.gz

  2. Need to modify:
    
    
    a) edit kma indexing and alignment so that index files would not need to be copped in data directory (currently not elegant)
 

## Install
wget for yaml

## Preprocessing assumptions
Fastq files have been filtered for quality, adaptor sequences, optical duplicates, and host reads

## Dependencies
R, Python, KMA, CSVTK, Kraken2 (version ?), pandas

## Inputs
- reads in paired fastq files (.gz extension only as of 2024.02.02)
- out directory

## Example Usage
main.nf --reads {dir}  --outdir {dir}

## Outputs
File Name            |     Explanation        
------------------|--------------
abundance_plot.png | bar plot of RPKM by gene
gene_abundance_table.tsv | summarized RPKM values for all samples
kma_all_joinned.tsv | contains read counts for RPKM calculations (numerator)
kraken_all_report.tsv | bacterial read count for all samples (denominator)

next up: add an explanation of each column of output file

## Formula
The formula to calculate gene relative abundances is given by:

\[ \frac{{\text{{Gene Reads}}}}{{\text{{Length of gene per kb}}}} \times (\text{{Bacteria Depth}}) \times 10^9 \]

This formula was adapted from Munk et al. 2022

## Citations
###Kraken
Wood, D.E., Lu, J. & Langmead, B. Improved metagenomic analysis with Kraken 2. Genome Biol 20, 257 (2019). [https://doi.org/doi-number](https://doi.org/10.1186/s13059-019-1891-0)
###Formula
Munk, P., Brinch, C., Møller, F.D. et al. Genomic analysis of sewage from 101 countries reveals global landscape of antimicrobial resistance. Nat Commun 13, 7251 (2022). [https://doi.org/doi-number](https://doi-org/10.1038/s41467-022-34312-7)

