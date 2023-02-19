#!/usr/bin/env bash
#SBATCH -p long
#SBATCH -J pyani
#SBATCH --cpus-per-task=64
#SBATCH --mem=64G
#SBATCH -N 1 
#SBATCH --ntasks-per-node=1

source miniconda3/bin/activate pyani

#Calculate ANI for a set of genome assemblies 

genome_dir=$1
outdir=$2

average_nucleotide_identity.py \
-i $genome_dir \
-o $outdir \
-m ANIm 


