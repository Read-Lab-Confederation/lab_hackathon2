# lab_hackathon2
Our project is to adapt Brooke's script for determining AMR reads per million from meta genomic data

Here is the [the google doc](https://docs.google.com/document/d/1a1NjFz8dDE0VPHwXtsbpe8BlrzB2dXce/edit) with project info.

[Lab meeting presentation](https://docs.google.com/presentation/d/1DkSnNEFyrNsgcvd66kASLn81-8fmhp7KXd5w81L1w_0/edit#slide=id.g17b97a35150_1_0) where we did git tutorial.

[Github site](https://github.com/Read-Lab-Confederation/github-collab-practice) for practising github pushes and pulls.

**Agenda for thurs 29th June 2023**

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

**Agenda for friday 29th June 2023**

1. Downsample the synthetic datasets?
2. One group work on Rscript
```Rscript abundanceTable.r --input all_kma_numerators_raw.tab --input_denominator all_kraken_bacteria.tab```
3. The other group work on converting the existing shell script to nextflow
4. Docker container?


**Andrei's work on the pipeline**
  1. Results and instructions:
    a) wrapped the whole pipeline in nextflow
    b) created new yaml for conda dependencies
    c) to install the environment type: conda env create -fhack2_nextflow.yaml
    d) to run the pipeline type:  bash -i runAll.sh
    e) current version uses local paths for references, so data/, bin/, main.nf, and runAll.sh should be in the same directory
    f) data/ directory can be downloaded from s3://transfer-files-emory/amrKma/data.tar.gz

  2. Need to modify:
    a) edit kma indexing and alignment so that index files would not need to be copped in data directory (currently not elegant)
 

