git clone https://github.com/Read-Lab-Confederation/lab_hackathon2.git
ls
cd lab_hackathon2/
ls
conda create -c bioconda -n hack2 nextflow kma kraken2 csvtk fastp -y
conda activate hack2
conda install conda-forge::r-base
conda install -y -c conda-forge r r-tidyr
conda install -y -c conda-forge r r-ggplot2
conda install conda-forge::matplotlib  
conda install conda-forge::pandas -y
mkdir data
data/
git clone https://github.com/tseemann/abricate.git
cd abricate/
mv db/ ../
cd ..
rm abricate/
ls
rm abricate/
cd ..
python db_merger.py
mv compiled_AMR_database.fasta data/
cd data
cp /mnt/tiramisu/emergent/projects/SEMAPHORE/data/fastqs/semaphore/microbiome/data_files/S.190905.00152_* ./
wget https://genome-idx.s3.amazonaws.com/kraken/k2_standard_08gb_20230605.tar.gz
tar -xvzf k2_standard_08gb_20230605.tar.gz
mkdir index
cd ..
nextflow run main.nf