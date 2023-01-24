#!/usr/bin/env bash
#SBATCH -p long
#SBATCH -J checkM
#SBATCH --mem=50G
#SBATCH --cpus-per-task=24

#Infiles should consist of a directory of genomes
#Outdir should be a location for results to be written

infile_dir=$1
outdir=$2

#Script to run checkM on genomes 

checkm lineage_wf \
-t 16 \
-x .fna \
$infile_dir \
$outdir


