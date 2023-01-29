#!/usr/bin/env bash
#SBATCH --mem=50G
#SBATCH -p long
#SBATCH -J pyani
#SBATCH --cpus-per-task=12

#Calculate ANI for a set of genome assemblies 

genome_dir=$1
outdir=$2

average_nucleotide_identity.py \
-i $genome_dir \
-o $outdir \
-m ANIm \
-g 

