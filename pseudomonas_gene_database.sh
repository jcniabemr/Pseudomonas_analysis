#!/usr/bin/env bash 
#SBATCH -p long 
#SBATCH -J PSGDB
#SBATCH --mem=4G
#SBATCH --cpus-per-task=4

####Script to document the commands used to create Ps gene database 

####Download DC3000 P. syringae pv. tomato DC3000 genome (NCBI accession no. NC_004578)
source miniconda3/bin/activate ncbi_datasets 
mkdir -p pseudomonas/gene_database/DC3000
cd pseudomonas/gene_database/DC3000
datasets download genome accession "GCF_000007805.1" --include genome,gff3,cds
unzip ncbi_dataset.zip
rm ncbi_dataset.zip

####Create NCBI database from DC3000 CDS